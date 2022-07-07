classdef FlowManager
    %FLOWMANAGER run on demand only the rquired parts of a recurrent algorithm's flow.
    %
    %Description:	
    %    The objective is to save the alrorithm runtime by running only parts that were
	%    not already calculated. This is done by sophisticated flowMng, and an internal
	%    cache. independent parts that might be calculated only once are stored in this 
	%    cache which is retreived on demand for the next parts. Each such part of the 
	%    flow is tagged by a "label", so it is called. The decision of which label to 
	%    run is made by preliminarily analyzing the length of all possible inputs for 
	%    each label. 
	%    
    %    The process is looked as sequential layers with dependent inputs and outputs.
    %    The entries of each layer is first registered. The entries are split into two 
    %    categories: those with one instance, and those with multiple instances. Since
    %    single instance entries require to be calculated only once. In the first run, 
    %    all the layers are run. All the entries are set.
    %
    %FLOWMANAGER methods:
    %    AddLabel(varargin) - Registers a label, and all its inputs.
    %    SetLabel(n) - Sets the number of the current label. 
    %    ok = start() - Checks if should run this part, or just retreiving the cache.
	%    SetReg(varargin) - Register the results of current label. 
	%    varargout = GetReg() - Retreives the results of current label. 
    %
    %Example:
	%    % The following algorithm calculates x=>y1=>y2=>y(n). But where y2 depends on n,
	%    % y1 is independent, thus might be calculated only once. Thus FLOWMANAGER calculates
	%    % it only in the first time, stores it in its cashe, and retreives it on the next
	%    % demands. 
    %    x = mydata;
    %    p1 = 2.5;
    %    p2 = 4.5;
    %    p3 = [12, 21, 345, 2];
    %    flowMng = FlowManager();
    %    flowMng = flowMng.AddLabel(p1);
    %    flowMng = flowMng.AddLabel(p2, p3);
    %  
    %    for n=1:numel(p3)	
    %        flowMng = flowMng.SetLabel(1);  % layer 1
    %        if (flowMng.start())  
    %            y1 = run_alg1(x, p1);  % runs some algorithm
    %            flowMng = flowMng.SetReg(y1);
    %        else
    %            y1 = flowMng.GetReg();
    %        end
    %  
    %        flowMng = flowMng.SetLabel(2);  % layer 2
    %        if (flowMng.start())
    %            y2 = run_alg2(y1, p2, p3(n));  % runs some algorithm
    %            flowMng = flowMng.SetReg(y2);
    %        else
    %            y2 = flowMng.GetReg();
    %        end
    %  
    %        y(n) = run_alg3(y2);  % layer 3 - runs some algorithm
    %    end
    %

    %    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
    %    Dept. Electrical & Computer Engineering, BGU Israel.
    %    This code is published under GNU GPLv3 license (see license in "LICENSE." file).
    
    properties
        enable
        reg
        N
        n
    end
    
    methods
        function obj = FlowManager()
            obj.N = 0;
            obj.n = 0;
        end
        
        function obj = AddLabel(obj, varargin)
            obj.N = obj.N + 1;
            obj.enable(obj.N) = any(cellfun(@numel, varargin)>1);
            obj.reg{obj.N} = [];
        end
        function obj = SetLabel(obj, n)
            obj.n = n;
        end
        function ok = start(obj)
            ok = any(arrayfun(@obj.startLabel, 1:obj.n));
        end
        function obj = SetReg(obj, varargin)
            obj.reg{obj.n} = varargin;
        end
        function varargout = GetReg(obj)
            varargout = obj.reg{obj.n};
        end       
    end
    methods (Access = private)
        function ok = startLabel(obj, n)
            ok = obj.enable(n) || isempty(obj.reg{n});
        end
    end
end

