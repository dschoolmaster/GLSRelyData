
#In this code, which is presented in Appendix S1, Box S1 we use the Inductive Causation Algorithm (IC-Algorithm) 
# given in Pearl (2009) section 2.5 titled "Recovering DAG Structures".
# This algorithm takes in a joint probability distribution of observed variables
# and uses the conditional independence relationships among them to
# determine a graph in which pairs of variables are concluded to be
# (1)  Non-adjacent (no edge between them, i.e. not connected by a link)
# (2)  Adjacent, but the direction of causation cannot be determined.
#      These pairs are left with a non-directed edge between them (i.e., a-b).
# (3)  Causally directed, where the pairs are connected a directed arrow 
#      (i.e., a->b) and are described 
#     by Pearl and Verma (1991, section 5) as "'genuine' causal relationships". 

#We apply this algorithm to a distribution consistent with that posited by
#  the SZC model for the relationships between Biodiversity,
#  Community composition and Trait composites. The goal is to see if 
#  the algorithm returns the causal relationships claimed by SZC, or 
#  if it finds that the relationships between community composition 
#  and the deterministic nodes Biodiversity and Trait composites is
#  non-causal, as claimed by GLS.

# To implement the algorithm, we make a distribution of observed variables 
# for three species (c1, c2, c3), Richness (R) and 
# Functional trait composites (i.e. trait diversity) (Q).

# For simplicity we treat:
# (1) The trait value of each species is treated as constant
# (2) We use 3 species to create variation in R, but do not show the results
#     for c3 to limit the size of the script. The tests and results are the same as with c1 and c2. 
# (3) We use familiar, linear, but approximate tests for the conditional independence
#     tests. Using more complicated, precise tests acts to strengthen the results.
# (4) We run the test 1000 times to demonstrate the stability of the distribution
#     to account for sampling errors.

#Set up variables to hold 1000 repeated simulations for our conditional
#independence tests
R_Q_null<-R_Q_C<-c1_R_null<-c2_R_null<-c1_Q_null<-c2_Q_null<-c1_c2_null<-rep(NA,1000)
c1_R_Z<-c2_R_Z<-c1_Q_Z<-c2_Q_Z<-matrix(NA,nrow = 1000,ncol=7)

for(i in 1:1000){
#Simulate 1000 communities
c1<-rnbinom(1000,mu = 10,size = .5)
c2<-rnbinom(1000,mu = 10,size = .5)
c3<-rnbinom(1000,mu = 10,size = .5)
#Set some trait value for each spp. (we treat these are constants for simplicity)
t1<-rnorm(1,10,3)
t2<-rnorm(1,10,3)
t3<-rnorm(1,10,3)

#Derive R and Q where Q is a weighted trait variation for each community

R<-apply(cbind(c1>0,c2>0,c3>0),1,sum)

#Create function to calculate Q which returns NA if c_x=0
Q_bits<-cbind(ifelse(c1==0,NA,c1*t1),ifelse(c2==0,NA,c2*t2),ifelse(c3==0,NA,c3*t3))
Q<-apply(Q_bits,1,function(x)sum(dist(x),na.rm = T))

#The resulting correlations are
example_cor_structure<-cor(cbind(c1,c2,c3,R,Q))

#We begin to apply the IC-Algorithm

#STEP 1 (Pearl 2009 Section 2.5):For each pair of variables test
#       if there is a set Z of variables
#       that makes them independent including the null set.
#       if they are non-independent given the null set and cannot be 
#       made independent by conditioning on any set Z then draw an un-directed
#       edge (a-b) between them.  

#Are R and Q independent by conditioning on Z={null}?
R_Q_null[i]<-summary(lm(R~Q))$coef["Q","Pr(>|t|)"]>0.05
#FALSE

#Can they be made independent by conditioning on composition (c1,c2,c3)?
R_Q_C[i]<-summary(lm(R~Q+(c1>0)+(c2>0)+(c3>0)))$coef["Q","Pr(>|t|)"]>0.05
#TRUE, so no edge is drawn between Q and R


#Are c1 and R independent by conditioning on Z={null}?
c1_R_null[i]<-summary(lm(c1>0~R))$coef["R","Pr(>|t|)"]>0.05
#FALSE

#Can c1 and R be made independent by conditioning on any set of vars?
c1_R_Z[i,]<-c(summary(lm(c1>0~R+Q+(c3>0)+(c2<0)))$coef["R","Pr(>|t|)"]>0.05,
              summary(lm(c1>0~R+Q+(c3>0)))$coef["R","Pr(>|t|)"]>0.05,
              summary(lm(c1>0~R+Q+(c2<0)))$coef["R","Pr(>|t|)"]>0.05,
              summary(lm(c1>0~R+(c3>0)+(c2<0)))$coef["R","Pr(>|t|)"]>0.05,
              summary(lm(c1>0~R+(c3>0)))$coef["R","Pr(>|t|)"]>0.05,
              summary(lm(c1>0~R+(c2<0)))$coef["R","Pr(>|t|)"]>0.05,
              summary(lm(c1>0~R+Q))$coef["R","Pr(>|t|)"]>0.05)
#FALSE so draw an un-directed edge between c1 and R (c1-R)

#Are c2 and R independent by conditioning on Z={null}?
c2_R_null[i]<-summary(lm(c2>0~R))$coef["R","Pr(>|t|)"]>0.05
#FALSE

#Can c2 and R be made independent by conditioning on any set of vars?
c2_R_Z[i,]<-c(summary(lm(c2>0~R+Q+(c3>0)+(c1>0)))$coef["R","Pr(>|t|)"]>0.05,
              summary(lm(c2>0~R+Q+(c3>0)))$coef["R","Pr(>|t|)"]>0.05,
              summary(lm(c2>0~R+Q+(c1>0)))$coef["R","Pr(>|t|)"]>0.05,
              summary(lm(c2>0~R+(c3>0)+(c1>0)))$coef["R","Pr(>|t|)"]>0.05,
              summary(lm(c2>0~R+(c3>0)))$coef["R","Pr(>|t|)"]>0.05,
              summary(lm(c2>0~R+(c1>0)))$coef["R","Pr(>|t|)"]>0.05,
              summary(lm(c2>0~R+Q))$coef["R","Pr(>|t|)"]>0.05)
#FALSE so draw an un-directed edge between c2 and R (c2-R)
    
#Are c1 and Q separated by Z={null}?
c1_Q_null[i]<-summary(lm(c1~Q))$coef["Q","Pr(>|t|)"]>0.05
#FALSE

#Can c1 and Q be made independent by conditioning on any set of vars?
c1_Q_Z[i,]<-c(summary(lm(c1~Q+c2+c3))$coef["Q","Pr(>|t|)"]>0.05,
              summary(lm(c1~Q+c2))$coef["Q","Pr(>|t|)"]>0.05,
              summary(lm(c1~Q+c3))$coef["Q","Pr(>|t|)"]>0.05,
              summary(lm(c1~Q+R))$coef["Q","Pr(>|t|)"]>0.05,
              summary(lm(c1~Q+c2+R))$coef["Q","Pr(>|t|)"]>0.05,
              summary(lm(c1~Q+R+c3))$coef["Q","Pr(>|t|)"]>0.05,
              summary(lm(c1~Q+c2+c3+R))$coef["Q","Pr(>|t|)"]>0.05)
#FALSE so draw an un-directed edge between c1 and R (c1-Q)

#Are c2 and Q separated by Z={null}?
c2_Q_null[i]<-summary(lm(c2~Q))$coef["Q","Pr(>|t|)"]>0.05

#Can c2 and Q be made independent by conditioning on any set of vars?
c2_Q_Z[i,]<-c(summary(lm(c2~Q+c1+c3))$coef["Q","Pr(>|t|)"]>0.05,
              summary(lm(c2~Q+c1))$coef["Q","Pr(>|t|)"]>0.05,
              summary(lm(c2~Q+c3))$coef["Q","Pr(>|t|)"]>0.05,
              summary(lm(c2~Q+R))$coef["Q","Pr(>|t|)"]>0.05,
              summary(lm(c2~Q+c1+R))$coef["Q","Pr(>|t|)"]>0.05,
              summary(lm(c2~Q+R+c3))$coef["Q","Pr(>|t|)"]>0.05,
              summary(lm(c2~Q+c1+c3+R))$coef["Q","Pr(>|t|)"]>0.05)
#FALSE so draw an un-directed edge between c2 and Q (c2-Q)

#Are c2 and c1 separated by Z={null}?
c1_c2_null[i]<-summary(lm(c1~c2))$coef["c2","Pr(>|t|)"]>0.05
#TRUE
}

#Look at an example correlation matrix of the variables.
example_cor_structure

#Reminder of step 1 procedure:  
#  for each pair of variables, we test if they are independent given a null set
#  (which asks if they are correlated in the joint distribution). If they 
#  are, we ask if there is any set of the other variables that makes
#  them independent by conditioning. If such a set exists, we do nothing.
#  If such a set does not exist, we connect the variables with an un-directed
#  edge (a-b). Figuring out the direction of these edges, if possible, happens
#  in the later steps.

#Percent of sims where R and Q are independent given null 
# (should be <<0.05 if associated with no conditioning set)
sum(R_Q_null)/1000

#Percent of sims where R and Q are independent given C
# (should near >=0.95 if R and Q are independent given C)
sum(R_Q_C)/1000

#Percent of sims where c1 or c2 are of R or Q independent given null
# (these should be <<0.05 if associated with no conditioning set)
sum(c1_R_null)/1000
sum(c2_R_null)/1000
sum(c1_Q_null)/1000
sum(c2_Q_null)/1000

#Percent of sims the c1 and c2 are independent given null
# (should near >=0.95 if not associated with no conditioning set)
sum(c1_c2_null)/1000

#Percent of sims where c1 or c2 are of R or Q independent given other variables
# (should be <=0.05 if {c1,R,Q} and {c2,R,Q} are not independent given any 
# conditioning set)
sum(apply(c1_R_Z,1,any))/1000
sum(apply(c2_R_Z,1,any))/1000
sum(apply(c1_Q_Z,1,any))/1000
sum(apply(c2_Q_Z,1,any))/1000


# Results of step 1:
# No edge is drawn between R and Q.
# No edge is drawn between c1 and c2.
# Un-directed edges are drawn between {(c1-R),(c2-R),(c1-Q),(c2-Q)}.

#STEP 2 (Pearl 2009 Section 2.5): For each pair of non-adjacent variables a and b (i.e. with no edge between them)
# but with a common neighbor c, check to see if the common neighbor c is in the 
# conditioning set that makes them independent
#  (a) if true do nothing 
#  (b) if false, then add arrowheads pointing at c from a and b (i.e., a->c<-b)

#From step 1, We have 2 pairs of non-adjacent variables (R,Q) and (c1,c2)

#Let's check R and Q. They are non-adjacent in our list from step 1
#  i.e. {(c1-R),(c2-R),(c1-Q),(c2-Q)}. They share common neighbors c1 and c2.
# However, c1 and c2 are both in the conditioning set that renders them 
# independent. Therefore, we do nothing (step 2 condition (a)).

#Let's check c1 and c2. They are non-adjacent in our list from step 1
# {(c1-R),(c2-R),(c1-Q),(c2-Q)}. They share common neighbors R and Q.
# Since neither R nor Q are not in the conditioning set that renders them
# independent (i.e. the null set), we draw the arrows c1->R<-c2 and c1->Q<-c2.
# (step 2 condition (b))

#There are additional steps to the algorithm (see: Pearl 2009 Section 2.5),
#but since all of our un-directed edges from step 1 are accounted for,
#we can stop.


#This demonstrates the IC-Algorithm has returned the casual relationships between 
#composition, richness and trait diversity posited by SZC


#according to Pearl and Verma (1991, Section 5):
 # "The result of this procedure is a substructure
#  called core(P^) in which every marked uni-directed arrow X->Y 
#  stands for the statement 'X has a causal influence on Y 
#  (in all minimal latent structures consistent with the data').
#  We call these relationships 'genuine' causal influences..."
# where the ... above refers to a figure in Pearl and Verma (1991, Figure 1a).

#Optional -- Draw resulting structure

#If the igraph package is not installed
# un-comment the following line:
#install.packages("igraph")

library(igraph)
#Declare nodes
nodes<-c("R","Q","c1","c2","c3")
#Define links resulting from the IC-Algorithm
links<-cbind("from"=c(rep("c1",2),rep("c2",2),rep("c3",2)),"to"=rep(c("R","Q"),3))
#Use these to create a igraph object
net<-graph_from_data_frame(links,vertices = nodes,directed = T)    
#Plot it to the screen
plot(net,vertex.size=30)

####Citations####
#Pearl, J. & T.S. Verma. 1991. A theory of inferred causation.
#  Proceedings of the Second International Conference on 
#  Principles of Knowledge Representation and Reasoning. 441–452.

# Pearl, J. 2009. Causality. Second edition. Cambridge University Press, Cambridge, UK

