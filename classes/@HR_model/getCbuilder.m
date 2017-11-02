function str = getCbuilder(object)
    str = sprintf('HR_model(%e,%e,%e,%e,%e)',object.b,object.I,object.mu,object.s,object.x_rest);
end