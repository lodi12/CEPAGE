function str = getCbuilder(object)
    str = sprintf('Danner_model(%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e)',object.C,object.gNaP,object.ENa,object.gl,object.El,object.VhalfM,object.km, object.Vhalfh,object.kh,object.tau0,object.tauMax,object.VhalfTau,object.kTau,object.Di, object.gSynE,object.EsynE);
end


