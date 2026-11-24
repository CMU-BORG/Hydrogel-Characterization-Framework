%% Log Transformation for strictly positive variables
%{
    This class implements the functions needed to allow for sampling of
    log-transformed variables in the SMC and ABC-SMC classes. This class
    utilizes a log transformation to convert strictly positive values to
    the range (-inf,inf).

    Last Updated: MJB 10/24/2026
%}

classdef LogitTransform
    %{


    %}

    properties
        a   % minimum value
        b   % maximum value
    end

    methods
        function obj = LogitTransform(a,b)
            obj.a = a;
            obj.b = b;
        end

        function logit_p = transform(obj,p)
            logit_p = log((p - obj.a)./(obj.b - p));
        end

        function p = inv_transform(obj,logit_p)
            p = obj.a + (obj.b - obj.a)./(1 + exp(-logit_p));
        end

        function J = jacobian(obj,p)
            J = (obj.b - obj.a)./( (p-obj.a).*(obj.b-p) );
        end

        function J_inv = jacobian_inv(obj, logit_p)
            J_inv = 1./obj.jacobian(obj.inv_transform(logit_p));
        end
    end
end