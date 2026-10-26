%% Log Transformation for strictly positive variables
%{
    This class implements the functions needed to allow for sampling of
    log-transformed variables in the SMC and ABC-SMC classes. This class
    utilizes a log transformation to convert strictly positive values to
    the range (-inf,inf).

    Last Updated: MJB 10/24/2026
%}

classdef LogTransform
    %{


    %}

    properties


    end

    methods
        function obj = LogTransform()

        end

        function log_p = transform(obj,p)
            log_p = log(p);
        end

        function p = inv_transform(obj,log_p)
            p = exp(log_p);
        end

        function J = jacobian(obj,p)
            J = abs(1./p);
        end

        function J_inv = jacobian_inv(obj, exp_p)
            J_inv = 1./obj.jacobian(obj.inv_transform(exp_p));
        end
    end
end