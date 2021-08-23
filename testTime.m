clc
clear
close all
f=figure('menubar','none');
hold on;
for i=1:60
    tmp=i*pi/30;
    plot([0,5*sin(tmp)],[0,5*cos(tmp)],'.k');
end
for i=1:12
    tmp=i*pi/6;
    plot([0,5*sin(tmp)],[0,5*cos(tmp)],'*r');
end
axis([-6,6,-6,6])
axis equal
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')


try
    while 1
        time=floor(clock);
        hour=time(4);
        min=time(5);
        sec=time(6);
        
        argHour=(hour+min/60)*pi/6;
        hHour=plot([0,3*sin(argHour)],[0,3*cos(argHour)],'r');
        
        argMin=(min+sec/60)*pi/30;
        hMin=plot([0,4*sin(argMin)],[0,4*cos(argMin)],'b');
        
        argSec=sec*pi/30;
        hSec=plot([0,4.5*sin(argSec)],[0,4.5*cos(argSec)],'k');
        
        set(f,'Name',['北京时间',num2str(hour),' : ',num2str(min),' : ',num2str(sec)])
        
        drawnow;
        pause(1);
        delete([hHour,hMin,hSec])
        
%          title('北京时间 ',num2str(hour),' : ',num2str(min))
    end
catch
    
end