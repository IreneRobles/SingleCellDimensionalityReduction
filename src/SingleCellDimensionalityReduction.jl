

using DimensionalityReduction
using SingleCellExperimentJulia
using PyPlot
using JLD


function DimensionalityReduction.dimensionality_reduction(method::String, 
        sceexp::SingleCellExperimentJulia.SingleCellExperiment, 
        outdim::Int; 
        assay = "CPM",
        k=100, 
        t::Int=1, 
        ɛ::Float64=1.0, 
        noise=0.001:0.00001:0.002, 
        initial_dims = size(sceexp.rownames)[1], 
        iterations=1000, 
        perplexity=20,
        varcomponents = 0
    )
    
    matrix = Array(sceexp.assays[assay]')
    matrix =   float(matrix)
    
    if method == "tSNE" && varcomponents > 0
        # Perform PCA noise reduction
          matrix = Array(DimensionalityReduction.dimensionality_reduction("PCA", matrix, varcomponents,
    k=k, t=t, ɛ=ɛ, noise=noise, initial_dims = initial_dims, iterations=iterations, perplexity=perplexity))
         
         red_dims = DimensionalityReduction.dimensionality_reduction(method, matrix, outdim,
    k=k, t=t, ɛ=ɛ, noise=noise, initial_dims = initial_dims, iterations=iterations, perplexity=perplexity)
        
        else
    
        red_dims = DimensionalityReduction.dimensionality_reduction(method, matrix, outdim,
    k=k, t=t, ɛ=ɛ, noise=noise, initial_dims = initial_dims, iterations=iterations, perplexity=perplexity)
        end
    
    # Store Variances associated to component if method is PCA
    if method == "PCA"
        pca =  DimensionalityReduction.dimensionality_reduction_model(method, matrix, outdim)
        variances = pca.prinvars ./ pca.tvar
        sceexp.metadata[string(method, "_assay=", assay,"_outdims=", outdim, "_VarianceComponents")] = variances
    end
    
        
    # Save the dimensionality reduction into the sce object
    if method == "tSNE"
        sceexp.reducedDimNames[string(method, "_assay=", assay,"_outdims=", outdim, "_perplexity=", perplexity, "_PCAcomps=", varcomponents)] = red_dims
    elseif method == "Isomap"
        sceexp.reducedDimNames[string(method, "_assay=", assay,"_outdims=", outdim, "_k=", k)] = red_dims
    else
        sceexp.reducedDimNames[string(method, "_assay=", assay,"_outdims=", outdim)] = red_dims 
    end
    return sceexp
    
end
    
        
    
        
        
######################
## PLOTTNG ##
#############

export plot_PCA_Variance
        
function plot_PCA_Variance(sceexp; ncomps = 15, assay = "CPM")
    sceexp = DimensionalityReduction.dimensionality_reduction("PCA", sceexp, ncomps, assay = assay)
    vars = sceexp.metadata[string("PCA", "_assay=", assay,"_outdims=", ncomps, "_VarianceComponents")]
    bar(1:length(vars), vars*100, 0.4, align="center")
    
    title(string("PCA $assay"))
    ylabel("% Variance")
    xlabel("Component")
end

