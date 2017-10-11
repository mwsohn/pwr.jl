module pwr

using Distributions, Roots, DataStructures

export cohenES, ESh, ESw1, ESw2
export powerAnovaTest, samplesizeAnovaTest, effectsizeAnovaTest, alphaAnovaTest
export powerChisqTest, samplesizeChisqTest, effectsizeChisqTest, alphaChisqTest
export powerF2Test, samplesizeF2Test, effectsizeF2Test, alphaF2Test
export powerNormTest, samplesizeNormTest, effectsizeNormTest, alphaNormTest
export powerPTest, samplesizePTest, effectsizePTest, alphaPTest
export powerTwoPTest, samplesizeTwoPTest, effectsizeTwoPTest, alphaTwoPTest
export powerTwoP2nTest, samplesizeTwoP2nTest, effectsizeTwoP2nTest, alphaTwoP2nTest
export powerRTest, samplesizeRTest, effectsizeRTest, alphaRTest
export powerT2nTest, samplesizeT2nTest, effectsizeT2nTest, alphaT2nTest
export powerTTest, samplesizeTTest, effectsizeTTest, alphaTTest
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
include("TTest.jl")
include("TwoP2nTest.jl")
include("TwoPTest.jl")
include("PowerPlot.jl")

end
