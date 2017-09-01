module pwr

using Distributions, Roots, DataStructures

export ESh, ESw1, ESw2
export powerAnovaTest, samplesizeAnovaTest, effectsizeAnovaTest, alphaAnovaTest
export powerChisqTest, samplesizeChisqTest, effectsizeChisqTest, alphaChisqTest
export powerPTest, samplesizePTest, effectsizePTest, alphaPTest
export powerT2Test, samplesizeT2Test, effectsizeT2Test, alphaT2Test
export powerTTest, samplesizeTTest, effectsizeTTest, alphaTTest
export power2p2nTest, samplesize2p2nTest, effectsize2p2nTest, alpha2p2nTest
export power2pTest, samplesize2pTest, effectsize2pTest, alpha2pTest

include("utils.jl")
include("ESh.jl")
include("ESw1.jl")
include("ESw2.jl")
include("AnovaTest.jl")
include("ChisqTest.jl")
include("PTest.jl")
include("T2Test.jl")
include("TTest.jl")
include("Twop2nTest.jl")
include("TwopTest.jl")

end
