function [xreset,object] = resetStates(object,t,x)

xreset(1) = object.c;
xreset(2) = x(2)+object.d;

end

