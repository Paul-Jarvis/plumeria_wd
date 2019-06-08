%octave script that plots output from plumeria fortran model

%*****************************************************************************
%THE FOLLOWING PARAMETERS NEED TO BE SET BEFORE RUNNING PLOTEMUP
% fprintf(1,'\n');
% fprintf(1,'Note: before using plotemup.m, you need to modify the first few\n');
% fprintf(1,'lines of the code to indicate the number and name(s) of the output file(s)\n');
% fprintf(1,'that will be read and plotted\n');
% fprintf(1,'Press <enter> to continue\n');
% pause;
n_outfiles = 1;                                      %number of files to read (<=2)
outputfile1 = 'output/Redoubt_200903231200Z_out.txt'
file1label  = 'Redoubt';                             %name of file 2 for plot legend
outputfile2 = 'output/5ms_east.txt';     %name & location of first file
file2label  = 'east';                             %name of file 2 for plot legend
outputfile3 = 'output/5ms_north.txt';     %name & location of first file
file3label  = 'north';                             %name of file 2 for plot legend
outputfile4 = 'output/5ms_NE.txt';     %name & location of first file
file4label  = 'NE';                             %name of file 2 for plot legend

%******************************************************************************
%WRITE OUT INPUT VALUES
fprintf(1,'number of output files to read = %d\n',n_outfiles);
fprintf(1,'output file 1=%s\n',outputfile1);
if n_outfiles==2,
   fprintf(1,'output file 2=%s\n',outputfile2);
end

%CLEAR THE OLD PLOT(S)
clear plot;
clear legend;

%******************************************************************************
%OPEN THE FIRST OUTPUT FILE
fprintf(1,'opening %s \n',outputfile1);
title=char(80);
fid=fopen(outputfile1);                          %open file

%Look for the title to the output table
while (length(title)<9) || ~strcmp(title(1:10),'**********') && ~feof(fid)
     title = fgetl(fid);  %read in a line, looking for the start of the table
     %fprintf(1,'%s\n',title);
     %pause
end

%If we didn't reach the end of the file, start reading numbers
if ~feof(fid)
    title=fgetl(fid);     %read the first line of the table titles
    title=fgetl(fid);     %read the second line of the table titles
    title=fgetl(fid);     %Now read the first line of numbers
    i=1;
    while ~strcmp(title(1:10),'**********')  && ~feof(fid)
    	numbers(:,i)=str2num(title);        %convert the line to numbers
        title=fgetl(fid);                   %read the next line
        i=i+1;
    end
end

    %Start assigning numbers to variables
%    if i==1,
       inum    = numbers(1,:);
       z       = numbers(3,:)/1000;
       x       = numbers(4,:)/1000;
       y       = numbers(5,:)/1000;
       u       = numbers(12,:);
       r       = numbers(13,:)/1000;
       T_mix   = numbers(14,:);
       T_air   = numbers(15,:);
       rho_mix = numbers(16,:);
       rho_air = numbers(17,:);
       time    = numbers(18,:)/60;
       p_air   = numbers(19,:);

       dist    = sqrt(x.^2 + y.^2);

fclose(fid);

maxz=max(z);
maxx=max(dist);
maxdist=max(dist);

%******************************************************************************
%OPEN THE SECOND OUTPUT FILE
if n_outfiles>=2,
   fprintf(1,'opening %s\n',outputfile2);
   fid2=fopen(outputfile2);    %open file
   fopen(fid2);
   title=char(80);
   %Look for the beginning of the output table
   while (length(title)<9) || ~strcmp(title(1:10),'**********') && ~feof(fid2)
        title = fgetl(fid2);  %read in a line, looking for the start of the table
   end

   %If we didn't reach the end of the file, start reading numbers
   if ~feof(fid2)
       title=fgetl(fid2);     %read the first line of the table titles
       title=fgetl(fid2);     %read the second line of the table titles
       title=fgetl(fid2);     %Now read the first line of numbers
       i=1;
       fprintf(1,'reading lines\n');
       while ~strcmp(title(1:10),'**********')  && ~feof(fid2)
           numbers2(:,i)=str2num(title);        %convert the line to numbers
           title=fgetl(fid2);                   %read the next line
           %fprintf(1,'%s\n',title);
           i=i+1;
       end
   end

   inum2   = numbers2(1,:);
   z2      = numbers2(3,:)/1000;
   x2      = numbers2(4,:)/1000;
   y2      = numbers2(5,:)/1000;
   u2      = numbers2(12,:);
   r2      = numbers2(13,:)/1000;
   T_mix2  = numbers2(14,:);
   T_air2  = numbers2(15,:);
   rho_mix2= numbers2(16,:);
   rho_air2= numbers2(17,:);
   time2   = numbers2(18,:)/60;
   p_air2  = numbers2(19,:);
   fclose(fid2);

   dist2   = sqrt(x2.^2 + y2.^2);
   
   if max(z2)>maxz,
       maxz=max(z2);
   end
   if max(dist2)>maxdist,
       maxdist=max(dist2);
   end
end

%******************************************************************************
%OPEN THE THIRD OUTPUT FILE
if n_outfiles>=3,
   fprintf(1,'opening %s\n',outputfile3);
   fid3=fopen(outputfile3);    %open file
   fopen(fid3);
   title=char(80);
   %Look for the beginning of the output table
   while (length(title)<9) || ~strcmp(title(1:10),'**********') && ~feof(fid3)
        title = fgetl(fid3);  %read in a line, looking for the start of the table
   end

   %If we didn't reach the end of the file, start reading numbers
   if ~feof(fid3)
       title=fgetl(fid3);     %read the first line of the table titles
       title=fgetl(fid3);     %read the second line of the table titles
       title=fgetl(fid3);     %Now read the first line of numbers
       i=1;
       fprintf(1,'reading lines\n');
       while ~strcmp(title(1:10),'**********')  && ~feof(fid3)
           numbers3(:,i)=str2num(title);        %convert the line to numbers
           title=fgetl(fid3);                   %read the next line
           %fprintf(1,'%s\n',title);
           i=i+1;
       end
   end

   inum3   = numbers3(1,:);
   z3      = numbers3(3,:)/1000;
   x3      = numbers3(4,:)/1000;
   y3      = numbers3(5,:)/1000;
   u3      = numbers3(12,:);
   r3      = numbers3(13,:)/1000;
   T_mix3  = numbers3(14,:);
   T_air3  = numbers3(15,:);
   rho_mix3= numbers3(16,:);
   rho_air3= numbers3(17,:);
   time3   = numbers3(18,:)/60;
   p_air3  = numbers3(19,:);
   fclose(fid3);

   dist3   = sqrt(x3.^2 + y3.^2);
   if max(z3)>maxz,
       maxz=max(z3);
   end
   if max(dist3)>maxdist,
       maxdist=max(dist3);
   end
end

%******************************************************************************
%OPEN THE FOURTH OUTPUT FILE
if n_outfiles>=4,
   fprintf(1,'opening %s\n',outputfile4);
   fid4=fopen(outputfile4);    %open file
   fopen(fid4);
   title=char(80);
   %Look for the beginning of the output table
   while (length(title)<9) || ~strcmp(title(1:10),'**********') && ~feof(fid4)
        title = fgetl(fid4);  %read in a line, looking for the start of the table
   end

   %If we didn't reach the end of the file, start reading numbers
   if ~feof(fid4)
       title=fgetl(fid4);     %read the first line of the table titles
       title=fgetl(fid4);     %read the second line of the table titles
       title=fgetl(fid4);     %Now read the first line of numbers
       i=1;
       fprintf(1,'reading lines\n');
       while ~strcmp(title(1:10),'**********')  && ~feof(fid4)
           numbers4(:,i)=str2num(title);        %convert the line to numbers
           title=fgetl(fid4);                   %read the next line
           %fprintf(1,'%s\n',title);
           i=i+1;
       end
   end

   inum4   = numbers4(1,:);
   z4      = numbers4(3,:)/1000;
   x4      = numbers4(4,:)/1000;
   y4      = numbers4(5,:)/1000;
   u4      = numbers4(12,:);
   r4      = numbers4(13,:)/1000;
   T_mix4  = numbers4(14,:);
   T_air4  = numbers4(15,:);
   rho_mix4= numbers4(16,:);
   rho_air4= numbers4(17,:);
   time4   = numbers4(18,:)/60;
   p_air4  = numbers4(19,:);
   fclose(fid4);

   dist4   = sqrt(x4.^2 + y4.^2);
   if max(z4)>maxz,
       maxz=max(z4);
   end
   if max(dist4)>maxdist,
       maxdist=max(dist4);
   end
end

maxdist=1.1*maxdist;
maxz   =1.1*maxz;

%******************************************************************************
%START PLOTTING
subplot(1,5,1)                                         %velocity
   plot(dist,z), hold on;
   if n_outfiles==1,
      legend(file1label);
    elseif n_outfiles==2,
      plot(dist2,z2,'r');
      legend(file1label,file2label)
    elseif n_outfiles==3,
      plot(dist2,z2,'r');
      plot(dist3,z3,'k');
      legend(file1label,file2label,file3label)
    elseif n_outfiles==4,
      plot(dist2,z2,'r');
      plot(dist3,z3,'k');
      plot(dist4,z4,'m');
      legend(file1label,file2label,file3label,file4label)
   end
   set(gca,'XTick',[-1:1:8]);
   xlabel('x, m/s');
   ylabel('z above vent, km');
   axis([-1 maxdist 0 maxz])
subplot(1,5,2)                                         %temperature
   plot(T_mix,z), hold on;
   if n_outfiles==2,
      plot(T_mix2,z2,'r');
    elseif n_outfiles==3,
      plot(T_mix2,z2,'r');
      plot(T_mix3,z3,'k');
    elseif n_outfiles==4,
      plot(T_mix2,z2,'r');
      plot(T_mix3,z3,'k');
      plot(T_mix4,z4,'m');
   end
   xlabel('T, K');
   axis([200 1200 0 maxz]);
   set(gca,'XTick',[200:400:1200],'YTickLabel',[]);
subplot(1,5,3)                                         %radius
   plot(r,z), hold on;
   if n_outfiles==2,
      plot(r2,z2,'r');
    elseif n_outfiles==3,
      plot(r2,z2,'r');
      plot(r3,z3,'k');
    elseif n_outfiles==4,
      plot(r2,z2,'r');
      plot(r3,z3,'k');
      plot(r4,z4,'m');
   end
   xlabel('r, km')
   axis([0 3 0 maxz]);
   set(gca,'YTickLabel',[]);
subplot(1,5,4)                                         %density
   plot(rho_mix,z), hold on;
   if n_outfiles==2,
      plot(rho_mix2,z2,'r');
    elseif n_outfiles==3,
      plot(rho_mix2,z2,'r');
      plot(rho_mix3,z3,'k');
    elseif n_outfiles==4,
      plot(rho_mix2,z2,'r');
      plot(rho_mix3,z3,'k');
      plot(rho_mix4,z4,'m');
   end
   xlabel('rho, kg/m^3');
   axis([0 6 0 maxz]);
   set(gca,'YTickLabel',[]);
subplot(1,5,5)                                        %time
   plot(time,z), hold on;
   if n_outfiles==2,
      plot(time2,z2,'r');
    elseif n_outfiles==3,
      plot(time2,z2,'r');
      plot(time3,z3,'k');
    elseif n_outfiles==4,
      plot(time2,z2,'r');
      plot(time3,z3,'k');
      plot(time4,z4,'m');
   end
   xlabel('t, minutes');
   axis([0 5 0 maxz]);
   set(gca,'XTick',[0:1:5],'YTickLabel',[]);
print('dataplot.ps','-deps')     %print the figure to a postscript file.
