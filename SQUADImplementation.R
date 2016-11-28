library(deSolve)
library(reshape2)
library(ggplot2)
#______________________________________________________________________________________________________________________#
#The following R code summarizes the implementation of the SQUAD method and presents several examples simulating 
#the effect of small and high intensity perturbations over the nodes of a regulatory network as described in 
#Fig. 4 of the Chapter: "The SQUAD method for the qualitative analysis of regulatory networks".
#______________________________________________________________________________________________________________________#

#__________________________________________NETWORK DEFINITION__________________________________________________________#
#First, define the regulatory network structure by postulating the corresponding logic rules and ODEs 
#for all nodes of the network. The SQUAD() function generates an standardized differential equation for each node. 
#The resulting ODEs system will then be numerically solved with the "ode" function of the "deSolve" package. 
#See "deSolve" package help for more information regarding ODEs system construction.
#______________________________________________________________________________________________________________________#

network<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    #Define the Fuzzy logic rules summarizing the regulatory information for each node
    #Classic logic operators AND, OR and NOT are translated into its Fuzzy logic equivalent by 
    #the use of the min(), max(), and 1-f() functions, respectively. 
    w_A<- min(max(X, A),1-B)
    w_B<- min(max(Y,B),1-A)
    w_X<- min(A,1-Z)
    w_Y<- min(B,1-A)
    w_Z<- B  
    #Model equations of the ODEs system using SQUAD
    dA<- SQUAD(A, w_A, gamma, h)
    dB<- SQUAD(B, w_B, gamma, h)
    dX<- SQUAD(X, w_X, gamma, h)
    dY<- SQUAD(Y, w_Y, gamma, h)
    dZ<- SQUAD(Z, w_Z, gamma, h)
    #ODE derivatives describing the dynamics of the nodes are returned as a list
    return(list(c(dA, dB, dX, dY, dZ)))})}


# Define SQUAD generic function to generate an ODE for each node of the regulatory network
SQUAD<-function(x,w,gamma,h){
  val<- ((-exp(0.5*h) + exp(-h*(w-0.5))) / ((1-exp(0.5*h)) * (1+exp(-h*(w-0.5))))) - (gamma*x)
  return(val)
}
#__________________________________GENERIC FUNCTIONS___________________________________________________________________#
#The perturbNodes() function creates an 'event' structure as implemented in the deSolve package to apply 
#the desired perturbations to the ODEs system at specific times. 
#Each perturbation may attain multiple values of strength and duration.
#See the 'events' section from the deSolve package 
#(https://cran.r-project.org/web/packages/deSolve/deSolve.pdf) for more information.
#______________________________________________________________________________________________________________________#

perturbNodes<-function(nodes, at.times, duration, intensity, time.step.size){
    pulses<-vars<-values<-NULL
      for(i in 1:length(nodes)){
        #Adjust the interval of each perturbation by its duration in time steps.
        interval<-duration[i]*time.step.size
        pulse<-seq(at.times[i], at.times[i] + interval, by=time.step.size)
        pulses<-c(pulses,pulse)
        vars<-c(vars,rep(nodes[i],length(pulse)))
        values<-c(values, rep(intensity[i], length(pulse)))
        }          
  return(data.frame(var=vars, time=pulses, value=values, method=rep("rep", length(pulses))))
}


#Time course plot function for the ODEs system simulation
plotPerturbation<-function(data, title, xlabel, ylabel){
  data<-melt(as.data.frame(data), id.vars="time")
  p<-ggplot(data = data, aes(x=time,y=value, colour=variable))
  p<-p+geom_point(size=1.0,alpha=0)+geom_path(size=1.5,alpha=1.0)
  p<- p + scale_colour_brewer(palette = "Set1")
  p<-p+theme_tufte()
  p<-p+scale_y_continuous(limits=c(0,1.0))+scale_x_continuous(limits=c(0,30))
  p<-p+labs(list(title=title,x=xlabel,y=ylabel, colour=""))
  return(p)}

#_______________________________________GENERAL SIMULATION PARAMETERS________________________________________#
# Define network parameters
parameters <- c(h = 50,gamma = 1)

#Time step size for the ODE system
time.step.size<-0.01

#Time interval to analyze the dynamic behavior of the regulatory network
time <- seq(0,30,by=time.step.size)

#Define the initial state from which to start the simulation
initial.state<-c(A=0, B=0, X=0, Y=0, Z=0)

#________________________DEFINE PERTURBATIONS WITH VARYING STRENGTH AND DURATION_______________________________________#
#   To generate a perturbation of the regulatory network use the perturbNodes() function:
#   Specify which 'nodes' of the network to perturb (as a character string), 
#   the times at which the perturbations will be applied ('at.times') which will be adjusted according     
#   to the time step size of the ODEs system,
#   The 'duration' (how many time steps) and the 'intensity' values for each perturbation.
#   EXAMPLES:
#   The following examples are based on Fig. 4 from the Chapter: "The SQUAD method for the qualitative analysis of regulatory networks".
#   EXAMPLE 1: Effect of low intensity perturbations. Simulate 'low' intensity perturbations with values 0< intensity < 0.5,  
#   for the "X" and "Y" nodes, respectively.   
#______________________________________________________________________________________________________________________#
perturbation_small<-perturbNodes(nodes=c("X", "Y"),at.times=c(10,20), duration=c(1,1), intensity=c(0.25,0.25), time.step.size)

#Simulate the dynamics of the regulatory network from an initial state and under the perturbation events previously generated with the preturbNodes() function.
#Dynamical simulation of the ODEs system using SQUAD
dynamics.small.perturbation <- ode(y=initial.state,times=time,func=network,parms=parameters, events=list(data=perturbation_small), rtol=10e-6, atol=10e-6)

#Plot the time course dynamics of the regulatory network for the example of small intensity perturbations:
colnames(dynamics.small.perturbation)<-c("time", "Node A", "Node B", "In X", "In Y", "Out Z")
plotPerturbation(dynamics.small.perturbation, title="Small perturbation example", xlabel = "time", ylabel = "Level of activation")

#_____EXAMPLE 2: Effect of high intensity perturbations with distinct order of appearance___________________________#
#Perturbations with different order od appearence and intensity wich may have an effect over the dynamical behavior 
#of regulatory networks can also be simulated with SQUAD. In the following example the input node "X" precedes a perturbation of the node "Y". Both signals attain "high" activation values i.e., they have a value of 1.0.
perturbation_X_first<-perturbNodes(nodes=c("X", "Y"),at.times=c(10,20), duration=c(1,1), intensity=c(1.0,1.0), time.step.size)
#Dynamical simulation of the ODEs system using SQUAD
dynamics.X.first  <- ode(y=initial.state,times=time,func=network,parms=parameters, events=list(data=perturbation_X_first), rtol=10e-6, atol=10e-6)
#Plot the effect of a high intensity perturbation of the node "X" preceding a perturbation of node "Y".
colnames(dynamics.X.first)<-c("time", "Node A", "Node B", "In X", "In Y", "Out Z")
plotPerturbation(dynamics.X.first, title="Perturbation example (X first)", xlabel = "time", ylabel = "Level of activation")

#Simulate the same perturbations as above but in different order of appearance, first a perturbation of node "Y" and then node "X". 
perturbation_Y_first<-perturbNodes(nodes=c("Y", "X"),at.times=c(10,20), duration=c(1,1), intensity=c(1.0,1.0), time.step.size)
#Dynamical simulation of the ODEs system using SQUAD
dynamics.Y.first  <- ode(y=initial.state,times=time,func=network,parms=parameters, events=list(data=perturbation_Y_first), rtol=10e-6, atol=10e-6)
#Plot the effect of a high intensity perturbation of the node "Y" preceding a perturbation of node "X".
colnames(dynamics.Y.first)<-c("time", "Node A", "Node B", "In X", "In Y", "Out Z")
plotPerturbation(dynamics.Y.first, title="Perturbation example (Y first)", xlabel = "time", ylabel = "Level of activation")
