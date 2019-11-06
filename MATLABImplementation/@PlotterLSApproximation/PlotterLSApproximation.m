classdef PlotterLSApproximation < handle

    properties
        m_LeastSquaresApproximations = [ClassicLeastSquares WeightedLeastSquares]';
    end
    
    properties (Hidden)
        handle_samplePoints
        handles_approximations
        handle_figure = figure
    end
    
	methods
        function obj = PlotterLSApproximation(varargin)
            for index = 1:nargin
                if strcmp(varargin{index}, 'figureHandle'); obj.handle_figure = varargin{index+1}; end
            end
        end
        
        plot_All(obj);
        show_SamplePoints(obj);
        hide_SamplePoints(obj);
        show_Approximation(obj, approximationIndices);
        hide_Approximation(obj, approximationIndices);
    end
    
end