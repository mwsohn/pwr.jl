using pwr, Plots

## Exercise 8.1 P. 357 from Cohen (1988)
#pwr.anova.test(f=0.28,k=4,n=20,sig.level=0.05)
t = pwr.AnovaTest(f=0.28,k=4,n=20,power=0.0)
plot(t)

## Exercise 6.1 p. 198 from Cohen (1988)
#pwr.2p.test(h=0.3,n=80,sig.level=0.05,alternative="greater")
t2 = pwr.TwoPTest(h=0.3,n=80,power=0.0,alternative="greater")
plot(t2)

## Exercise 7.3 p. 251
#pwr.chisq.test(w=0.346,df=(2-1)*(3-1),N=140,sig.level=0.01)
t3 = pwr.ChisqTest(w=0.346,df=(2-1)*(3-1),N=140,alpha=0.01,power=0.0)
plot(t3)

## Exercise 6.5 p. 203 from Cohen (1988)
#pwr.p.test(h=0.2,n=60,sig.level=0.05,alternative="two.sided")
t4 = pwr.PTest(h=0.2,n=60,alpha=0.05,alternative="two",power=0.0)
plot(t4,backend=plotly)


## medium effect size for the correlation test
#cohen.ES(test="r", size="medium")
cohenES(test="r", size="medium")

## sample size for a medium size effect in the two-sided correlation test
## using the conventional power of 0.80
#pwr.r.test(r=cohen.ES(test="r",size="medium")$effect.size,power=0.80, sig.level=0.05, alternative="two.sided")
t5 = pwr.RTest(n=0,r=cohenES(test="r",size="medium"),power=0.80, alpha=0.05, alternative="two.sided")
plot(t5)

## Exercise 6.5 p. 203 from Cohen
#h<-ES.h(0.5,0.4)
#h
h = ESh(0.5,0.4)

#pwr.p.test(h=h,n=60,sig.level=0.05,alternative="two.sided")
t6 = pwr.PTest(h=h,n=60,alpha=0.05,power=0.0,alternative="two")
plot(t6)


## Exercise 7.1 p. 249 from Cohen
#P0<-rep(1/4,4)
P0 = fill(1/4,4)

#P1<-c(0.375,rep((1-0.375)/3,3))
P1 = vcat(0.375,fill((1-0.375)/3,3))

#ES.w1(P0,P1)
w = ESw1(P0,P1)

#pwr.chisq.test(w=ES.w1(P0,P1),N=100,df=(4-1))
t7 = pwr.ChisqTest(w=ESw1(P0,P1),N=100,df=(4-1),power=0.0)
plot(t7)

# ES.w2
#prob<-matrix(c(0.225,0.125,0.125,0.125,0.16,0.16,0.04,0.04),nrow=2,byrow=TRUE)
prob = reshape([0.225,0.125,0.125,0.125,0.16,0.16,0.04,0.04],4,2)'

#ES.w2(prob)
ESw2(prob)

#pwr.chisq.test(w=ES.w2(prob),df=(2-1)*(4-1),N=200)
t8 = pwr.ChisqTest(w=ESw2(prob),df=(2-1)*(4-1),N=200,power=0.0)
plot(t8)


## Exercise 6.1 p. 198 from Cohen (1988)
#pwr.2p.test(h=0.3,n=80,sig.level=0.05,alternative="greater")
t9 = pwr.TwoPTest(h=0.3,n=80,alpha=0.05,alternative="greater",power=0.0)
plot(t9)

## Exercise 6.3 P. 200 from Cohen (1988)
#pwr.2p2n.test(h=0.30,n1=80,n2=245,sig.level=0.05,alternative="greater")
t10 = pwr.TwoP2nTest(h=0.30,n1=80,n2=245,alpha=0.05,power=0.0,alternative="greater")
plot(t10)

## Exercise 6.7 p. 207 from Cohen (1988)
#pwr.2p2n.test(h=0.20,n1=1600,power=0.9,sig.level=0.01,alternative="two.sided")
t11 =pwr.TwoP2nTest(h=0.20,n1=1600,power=0.9,alpha=0.01,alternative="two.sided")
plot(t11)

## Exercise 8.1 P. 357 from Cohen (1988)
#pwr.anova.test(f=0.28,k=4,n=20,sig.level=0.05)
t12 = pwr.AnovaTest(f=0.28,k=4,n=20,alpha=0.05)
plot(t12)

## Exercise 8.10 p. 391
#pwr.anova.test(f=0.28,k=4,power=0.80,sig.level=0.05)
t13 = pwr.AnovaTest(f=0.28,k=4,power=0.80,alpha=0.05)
plot(t13)

## Exercise 7.1 P. 249 from Cohen (1988)
#pwr.chisq.test(w=0.289,df=(4-1),N=100,sig.level=0.05)
t14 = pwr.ChisqTest(w=0.289,df=(4-1),N=100,alpha=0.05)
plot(t14)

## Exercise 7.3 p. 251
#pwr.chisq.test(w=0.346,df=(2-1)*(3-1),N=140,sig.level=0.01)
t15 = pwr.ChisqTest(w=0.346,df=(2-1)*(3-1),N=140,alpha=0.01)
plot(t15)

## Exercise 7.8 p. 270
#pwr.chisq.test(w=0.1,df=(5-1)*(6-1),power=0.80,sig.level=0.05)
t16 = pwr.ChisqTest(w=0.1,df=(5-1)*(6-1),power=0.80,alpha=0.05)
plot(t16)

## Exercise 9.1 P. 424 from Cohen (1988)
#pwr.f2.test(u=5,v=89,f2=0.1/(1-0.1),sig.level=0.05)
t17 = pwr.F2Test(u=5,v=89,f2=0.1/(1-0.1),alpha=0.05)

## Power at mu=105 for H0:mu=100 vs. H1:mu>100 (sigma=15) 20 obs. (alpha=0.05)
# sigma<-15
σ = 15
# c<-100
c = 100
# mu<-105
μ = 105
# d<-(mu-c)/sigma
d = (μ - c)/σ
# pwr.norm.test(d=d,n=20,sig.level=0.05,alternative="greater")
pwr.NormTest(d=d,n=20,alpha=0.05,alternative="greater")

## Sample size of the test for power=0.80
#pwr.norm.test(d=d,power=0.8,sig.level=0.05,alternative="greater")
pwr.NormTest(d=d,power=0.8,alpha=0.05,alternative="greater")

## Power function of the same test
#mu<-seq(95,125,l=100)
μ = collect(linspace(95,125,100))

#d<-(mu-c)/sigma
d = (μ - c)/σ

# plot(d,pwr.norm.test(d=d,n=20,sig.level=0.05,alternative="greater")$power,
# type="l",ylim=c(0,1))
# abline(h=0.05)
# abline(h=0.80)
power = [powerNormTest(d=x,n=20,alternative="greater") for x in d]
plot(d,power,ylim=[0,1],legend=false)
hline!([0.05,0.8])

## Power function for the two-sided alternative
# plot(d,pwr.norm.test(d=d,n=20,sig.level=0.05,alternative="two.sided")$power,
#       type="l",ylim=c(0,1))
# abline(h=0.05)
# abline(h=0.80)
plot(d,[powerNormTest(d=x,n=20,alpha=0.05,alternative="two.sided") for x in d],ylim=[0,1],legend=false)
hline!([0.05,0.8])

## Exercise 6.5 p. 203 from Cohen
# h<-ES.h(0.5,0.4)
h = ESh(0.5,0.4)
# h
# pwr.p.test(h=h,n=60,sig.level=0.05,alternative="two.sided")
pwr.PTest(h=h,n=60,alpha=0.05,alternative="two.sided")

## Exercise 6.8 p. 208
#pwr.p.test(h=0.2,power=0.95,sig.level=0.05,alternative="two.sided")
pwr.PTest(h=0.2,power=0.95,alpha=0.05,alternative="two.sided")

## Exercise 3.1 p. 96 from Cohen (1988)
#pwr.r.test(r=0.3,n=50,sig.level=0.05,alternative="two.sided")
pwr.RTest(r=0.3,n=50,alpha=0.05,alternative="two.sided")

#pwr.r.test(r=0.3,n=50,sig.level=0.05,alternative="greater")
pwr.RTest(r=0.3,n=50,alpha=0.05,alternative="greater")

## Exercise 3.4 p. 208
#pwr.r.test(r=0.3,power=0.80,sig.level=0.05,alternative="two.sided")
pwr.RTest(r=0.3,power=0.80,alpha=0.05,alternative="two.sided")

#pwr.r.test(r=0.5,power=0.80,sig.level=0.05,alternative="two.sided")
pwr.RTest(r=0.5,power=0.80,alpha=0.05,alternative="two.sided")

#pwr.r.test(r=0.1,power=0.80,sig.level=0.05,alternative="two.sided")
pwr.RTest(r=0.1,power=0.80,alpha=0.05,alternative="two.sided")



## One sample (power)
## Exercise 2.5 p. 47 from Cohen (1988)
#pwr.t.test(d=0.2,n=60,sig.level=0.10,type="one.sample",alternative="two.sided")
pwr.TTest(d=0.2,n=60,alpha=0.10,power=0.,sampletype="one.sample",alternative="two.sided")

## Paired samples (power)
## Exercise p. 50 from Cohen (1988)
#d<-8/(16*sqrt(2*(1-0.6)))
d = 8/(16*sqrt(2*(1-0.6)))
#pwr.t.test(d=d,n=40,sig.level=0.05,type="paired",alternative="two.sided")
pwr.TTest(d=d,n=40,alpha=0.05,sampletype="paired",power=0.,alternative="two.sided")

## Two independent samples (power)
## Exercise 2.1 p. 40 from Cohen (1988)
#d<-2/2.8
d = 2/2.8
#pwr.t.test(d=d,n=30,sig.level=0.05,type="two.sample",alternative="two.sided")
pwr.TTest(d=d,n=30,alpha=0.05,sampletype="two.sample",power=0.,alternative="two.sided")

## Two independent samples (sample size)
## Exercise 2.10 p. 59
#pwr.t.test(d=0.3,power=0.75,sig.level=0.05,type="two.sample",alternative="greater")
pwr.TTest(d=0.3,power=0.75,alpha=0.05,sampletype="two.sample",alternative="greater")

## Exercise 2.3 p. 437 from Cohen (1988)
#pwr.t2n.test(d=0.6,n1=90,n2=60,alternative="greater")
pwr.T2nTest(d=0.6,n1=90,n2=60,alternative="greater")
