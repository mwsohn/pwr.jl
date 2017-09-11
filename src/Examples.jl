using pwr

## Exercise 8.1 P. 357 from Cohen (1988)
#pwr.anova.test(f=0.28,k=4,n=20,sig.level=0.05)
pwr.AnovaTest(f=0.28,k=4,n=20,power=0.0)

## Exercise 6.1 p. 198 from Cohen (1988)
#pwr.2p.test(h=0.3,n=80,sig.level=0.05,alternative="greater")
pwr.TwopTest(h=0.3,n=80,power=0.0,alternative="greater")

## Exercise 7.3 p. 251
#pwr.chisq.test(w=0.346,df=(2-1)*(3-1),N=140,sig.level=0.01)
pwr.ChisqTest(w=0.346,df=(2-1)*(3-1),N=140,alpha=0.01,power=0.0)

## Exercise 6.5 p. 203 from Cohen (1988)
#pwr.p.test(h=0.2,n=60,sig.level=0.05,alternative="two.sided")
pwr.PTest(h=0.2,n=60,alpha=0.05,alternative="two",power=0.0)


## medium effect size for the correlation test
#cohen.ES(test="r", size="medium")
cohenES(test="r", size="medium")

## sample size for a medium size effect in the two-sided correlation test
## using the conventional power of 0.80
#pwr.r.test(r=cohen.ES(test="r",size="medium")$effect.size,power=0.80, sig.level=0.05, alternative="two.sided")
pwr.RTest(n=0,r=cohenES(test="r",size="medium").d["effectsize"],power=0.80, alpha=0.05, alternative="two.sided")

## Exercise 6.5 p. 203 from Cohen
#h<-ES.h(0.5,0.4)
#h
h = ESh(0.5,0.4)

#pwr.p.test(h=h,n=60,sig.level=0.05,alternative="two.sided")
pwr.PTest(h=h,n=60,alpha=0.05,power=0.0,alternative="two")


## Exercise 7.1 p. 249 from Cohen
#P0<-rep(1/4,4)
P0 = fill(1/4,4)

#P1<-c(0.375,rep((1-0.375)/3,3))
P1 = vcat(0.375,fill((1-0.375)/3,3))

#ES.w1(P0,P1)
w = ESw1(P0,P1)

#pwr.chisq.test(w=ES.w1(P0,P1),N=100,df=(4-1))
pwr.ChisqTest(w=ESw1(P0,P1),N=100,df=(4-1),power=0.0)

# ES.w2
#prob<-matrix(c(0.225,0.125,0.125,0.125,0.16,0.16,0.04,0.04),nrow=2,byrow=TRUE)
prob = reshape([0.225,0.125,0.125,0.125,0.16,0.16,0.04,0.04],4,2)'

#ES.w2(prob)
ESw2(prob)

#pwr.chisq.test(w=ES.w2(prob),df=(2-1)*(4-1),N=200)
pwr.ChisqTest(w=ESw2(prob),df=(2-1)*(4-1),N=200,power=0.0)


## Exercise 6.1 p. 198 from Cohen (1988)
#pwr.2p.test(h=0.3,n=80,sig.level=0.05,alternative="greater")
pwr.TwopTest(h=0.3,n=80,alpha=0.05,alternative="greater",power=0.0)
