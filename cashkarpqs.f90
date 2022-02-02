
      SUBROUTINE cashkarpqs(y,dydx,x,h_in,yscal,h_out,hnext,derivs)

!     Subroutine that integrates properties with elevation using the Cash-Karp
!     method:
!     Cash, J.R. and Karp, A.H., A variable order Runge-Kutta method for initial
!         value problems with rapidly varying right-hand sides, ACM Transactions
!         on Mathematical Software, v. 16, no. 3, pp. 201-222.

      use precis_param
      real(kind=ip) h_out,hnext,h_in,x,dydx(13),y(13),yerr(13), &
                           yscal(13), ytemp(13)
      !if stepsize is smaller than the number below, stop calculation
      real(kind=ip), PARAMETER   :: EPS_TINY  = 1.0e-17_ip     
      !if yerr/yscal exceeds the number below, reduce step size.
      real(kind=ip), PARAMETER   :: EPS_SMALL = 5.0e-4_ip            
      character*1 answer
      external derivs
      real(kind=ip) alpha, max_error,h,xnew
      real(kind=ip) maxerrnow(13)

      h=h_in
1     call cashkarp(y,dydx,x,h,ytemp,yerr,derivs)
      max_error=0.
      maxerrnow = 0._ip
      maxerrnow = abs(yerr/yscal)
      !write(6,*) 'maxerrnow=',maxerrnow
      !write(6,*) 'continue?'
      !read(5,'(a1)') answer
      !if (answer.eq.'n') stop
      max_error=maxval(abs(yerr/yscal))        !find maximum error, scaled to scale factors
      max_error=max_error/EPS_SMALL
      if(max_error.gt.1.)then
        alpha=max_error**0.25                  !adjustment following Dormund, J.R.
        h=0.9*h/min(alpha,10._ip)                  !"Numerical Methods for Differential Equations", p. 85
        xnew=x+h
        if((abs((xnew-x))/x).lt.EPS_TINY) then
            write(6,*) 'stepsize is approximately zero.  Continue (y/n)?'
            read(5,'(a1)') answer
            if (answer.eq.'n') stop
        end if
        goto 1
      else
        !increase h, but not by more than 5x
        hnext=min(0.9*h/(max_error**0.2),5.*h)
        h_out=h
        x=x+h
        y=ytemp
        return
      endif
      end subroutine cashkarpqs

!******************************************************************************

      SUBROUTINE cashkarp(y_in,dydx,x,h,yout,yerr,derivs)

!     Cash-carp ode solver

!     Uses the cash-karp method of solving ode's using the fifth-order Runge-Kutta solution
!     with embedded fourth-order solution:

!     y_out(i) = y_in(i) + c(1)* dydx(i) + c(2)*dydx2(i) + c(3)*dydx3(i) 
!                        + c(4)*dydx4(i) + c(5)*dydx5(i) + c(6)*dydx6(i)
!     where c(i) is a constant with values given in Cash and Karp, and dydx, dydx2 etc.
!           are the slopes of y with x evaluated at the following intermediate values
!           of x: x+h/5, x+3h/10, x+3h/5, x+h, 1+7h/8.

!     This is compared with the fourth-order Runge-Kutta solution:
!     ystar_out(i) = y_in(i) + cstar(1)* dydx(i) + cstar(2)*dydx2(i) + cstar(3)*dydx3(i)
!                            + cstar(4)*dydx4(i) + cstar(5)*dydx5(i) + cstar(6)*dydx6(i)
!     where cstar(i) are coefficients with different values.

!    The differences between the fourth-order and fifth-order solutions are used in the
!    quality step adjustment in cashkarpqs

      use precis_param

      real(kind=ip) h,x                                           !step size, initial x value
      real(kind=ip) dydx(13),dydx2(13),dydx3(13),dydx4(13), &
                                     dydx5(13),dydx6(13)            !dydx values at intermediate x
      real(kind=ip) y_in(13), yerr(13),yout(13), &
                                ytemp(13) , ystar_out(13)           !in, out & error values of y
      external derivs
      !constants
      real(kind=ip) b_21,b_31,b_32,b_41,b_42,b_43,b_51,b_52,b_53, &
      b_54,b_61,b_62,b_63,b_64,b_65,c_1,c_3,c_4,c_6,dc_1,dc_3,dc_4,dc_5,dc_6

      real(kind=ip) cstar_1,cstar_3,cstar_4,cstar_5,cstar_6
      character(len=1)  :: answer

      !Constants taken from p. 206 of Cash and Karp (ACM Transactions of on Mathematical
      !Software, v. 16, pp. 201-222)
      b_21=.2
      b_31=3./40.
      b_32=9./40.
      b_41=.3
      b_42=-.9
      b_43=1.2
      b_51=-11./54.
      b_52=2.5
      b_53=-70./27.
      b_54=35./27.
      b_61=1631./55296.
      b_62=175./512.
      b_63=575./13824.
      b_64=44275./110592.
      b_65=253./4096.
      c_1=37./378.
      c_3=250./621.
      c_4=125./594.
      c_6=512./1771.
      cstar_1=2825./27648.
      cstar_3=18575./48384.
      cstar_4=13525./55296.
      cstar_5=-277./14336.
      cstar_6=1./4.
      dc_1=c_1-2825./27648.
      dc_3=c_3-18575./48384.
      dc_4=c_4-13525./55296.
      dc_5=277./14336.
      dc_6=c_6-.25

      !Calculate y and dydx at intermediate values of x
      ytemp=y_in+b_21*h*dydx
      call derivs(x+h/5.,ytemp,dydx2)
      ytemp=y_in+h*(b_31*dydx+b_32*dydx2)                                    !at h/5
      call derivs(x+3.*h/10.,ytemp,dydx3)
      ytemp=y_in+h*(b_41*dydx+b_42*dydx2+b_43*dydx3)                         !at 3h/10
      call derivs(x+6.*h/10.,ytemp,dydx4)
      ytemp=y_in+h*(b_51*dydx+b_52*dydx2+b_53*dydx3+b_54*dydx4)              !at 6h/10
      call derivs(x+h,ytemp,dydx5)
      ytemp=y_in+h*(b_61*dydx+b_62*dydx2+b_63*dydx3+b_64*dydx4+b_65*dydx5)   !at h
      call derivs(x+h,ytemp,dydx5)
      ytemp=y_in+h*(b_61*dydx+b_62*dydx2+b_63*dydx3+b_64*dydx4+b_65*dydx5)   !at h
      call derivs(x+0.875*h,ytemp,dydx6)
      
      !Fifth order Runge-Kutta formula
      yout=y_in+h*(c_1*dydx+c_3*dydx3+c_4*dydx4+c_6*dydx6)
      !Embedded fourth-order Runge-Kutta formula
      ystar_out = y_in+h*(cstar_1*dydx+cstar_3*dydx3+cstar_4*dydx4+ cstar_5*dydx5+cstar_6*dydx6)
      
      !Difference between fourth-order and fifth-order solutions
      yerr = yout-ystar_out

      return
      
      end subroutine cashkarp

