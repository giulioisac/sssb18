function[distance]=dist_particles(x1,x2,length,Gconst,dx)
    distance=min(abs(x1-x2),length-abs(x1-x2));
    %if(distance>length -Gconst*dx)
    %    distance=length +dx- distance;
    %end
end