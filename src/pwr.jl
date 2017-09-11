module pwr

using Distributions, Roots, DataStructures

export cohenES, ESh, ESw1, ESw2
export powerAnovaTest, samplesizeAnovaTest, effectsizeAnovaTest, alphaAnovaTest
export powerChisqTest, samplesizeChisqTest, effectsizeChisqTest, alphaChisqTest
export powerF2Test, samplesizeF2Test, effectsizeF2Test, alphaF2Test
export powerNormTest, samplesizeNormTest, effectsizeNormTest, alphaNormTest
export powerPTest, samplesizePTest, effectsizePTest, alphaPTest
export powerRTest, samplesizeRTest, effectsizeRTest, alphaRTest
export powerT2nTest, samplesizeT2nTest, effectsizeT2nTest, alphaT2nTest
export powerT2Test, samplesizeT2Test, effectsizeT2Test, alphaT2Test
export powerTTest, samplesizeTTest, effectsizeTTest, alphaTTest
export power2p2nTest, samplesize2p2nTest, effectsize2p2nTest, alpha2p2nTest
export power2pTest, samplesize2pTest, effectsize2pTest, alpha2pTest
export plot

include("utils.jl")
include("cohenES.jl")
include("ESh.jl")
include("ESw1.jl")
include("ESw2.jl")
include("AnovaTest.jl")
include("ChisqTest.jl")
include("F2Test.jl")
include("NormTest.jl")
include("PTest.jl")
include("RTest.jl")
include("T2nTest.jl")
include("T2Test.jl")
include("TTest.jl")
include("Twop2nTest.jl")
include("TwopTest.jl")
include("PowerPlot.jl")

end
