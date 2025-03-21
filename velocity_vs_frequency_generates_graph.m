% First, load form the server the file velocity_summary or analog, 
% that contains C for OC, OA and Stroke
for i = 1:length(C)
M = C{i};
fc = M(:,2);
v = M(:,3);
[B,I] = sort(fc);
fc = B;
v = v(I);
if M(1,1) == 1% patients
plot(fc,v,'r');hold on;
else
plot(fc,v,'b');hold on;
end
end
xlabel('fc (Hz)');ylabel('v (m/s)');
title('red = patients; blue = subjects');
