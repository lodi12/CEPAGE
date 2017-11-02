function str = getCbuilder(object)
    str = sprintf('IZ_model(%e,%e,%e,%e,%e,%e,%e)',...
    object.a,object.b,object.c,object.d,object.I,object.gL,object.EL);
end