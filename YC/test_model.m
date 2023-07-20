function R = test_model(x,z,I)

[~,dummy] = system('rm output.out');

t = x(68:69);
tbot = x(67);
x = x(1:66);
x = x(:);

% convert interval percentages into soil moistures
x(53) = x(53)*x(52); %wangxi change
x(54) = x(51) + x(54)*(x(53)-x(51));

% make sure that the canopy top is not lower than the canopy bottom wangxi change
HVT = load('hvt.txt');
x(4) =x(4)*HVT; %wangxi change
x(5) = x(5) * x(4);%wangxi change

% save parameters
save('parms.txt','x','-ascii');

% multipli LAI and SAI
time_parms = load('time_parms_original.txt');
time_parms(1,:) = time_parms(1,:).*t(1);
time_parms(2,:) = time_parms(2,:).*t(2);
save('time_parms.txt','time_parms','-ascii');

% bottom temperature
save('tbot.txt','tbot','-ascii');

% run model
system('main.exe');

try

 fname = strcat('output.out');
 y = load(fname);
 y = y(:,[4,5,6,7,10]);

 for d = 1:5
  if length(I{d}>5000)
   II = I{d}; 
   R(d) = (y(II,d)-z(II,d))'*(y(II,d)-z(II,d));
   R(d) = sqrt(R(d)/length(II));
  else
   R(d) = 0/0;
  end
 end 

catch
 fprintf(' Failed ');
 for d = 1:5
  R(d) = 0/0;
 end
end


