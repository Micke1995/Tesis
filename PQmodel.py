import numpy as np
import matplotlib.pyplot as plt
from random import uniform,random,randint
from numpy import pi,exp,sin

class PQ:
    def __init__(self,Magnitud=1,Frecuency=60,Cicles=10,FS=16000,Amin=0.1,Amax=0.9,Phmin=-pi,Phmax=pi,PeriodoDisturbio=None,FinalDisturbio=0,InicioDisturbio=None):
    #### Parameteres of the Model   
        self.A = Magnitud
        self.f = Frecuency
        self.n = Cicles
        self.fs = FS
        self.PointsPerSignal = int((self.fs/self.f)*self.n)
        self.t=np.arange(0,(1/self.f)*self.n-(1/self.fs),1/self.fs)
        self.t1t2 = PeriodoDisturbio 
        self.PeriodLimit = FinalDisturbio
        self.Inicio = InicioDisturbio

    ### Random Parameters
        #General
        self.phase_min = Phmin
        self.phase_max = Phmax
        #Oscillatory transient and harmonics related
        self.theta_min = -pi
        self.theta_max = pi        
        #Sag,Swell,interrumption related
        self.Pmax = self.n-1
        self.Pmin = 0
        self.alpha_min = Amin
        self.alpha_max = Amax
        self.beta_min =0.1
        self.beta_max = 0.8
        self.rho_min=0.9
        self.rho_max=1
        #Transient related
        self.taPeriod_min = 1
        self.taPeriod_max = self.n-1
        self.psi_min = 0.222
        self.psi_max = 1.11
        self.Onems=round(0.001/(1/self.fs))# Number of points equivalent to 1ms
        # Flicker related
        self.ff_min = 8
        self.ff_max = 25
        self.lambda_min = 0.05
        self.lambda_max = 0.1
        # Oscillatory transient related
        self.tau_min = 0.008
        self.tau_max = 0.04
        self.fn_min = 300
        self.fn_max = 900
        self.periodMaxOT = self.n/3.33
        self.periodMinOT = 0.5
        self.pointsfithpartperiod=round(self.fs/(5*self.f))
        # Harmonics related
        self.alpha1=1
        self.alpha3_min=0.05
        self.alpha3_max=0.15
        self.alpha5_min=0.05
        self.alpha5_max=0.15
        self.alpha7_min=0.05
        self.alpha7_max=0.15
        # Notch related
        self.k_min=0.1
        self.k_max=0.4
        self.c=[1,2,4,6]# number of notches per cycle (possible values: 1,2,4 and 6)
        self.td_min=0
        self.tc_min=0
        self.tdminustc_min=0.01*(1/self.f)
        self.tdminustc_max=0.05*(1/self.f)

        
    
    def rphi(self):
        return uniform(self.phase_min,self.phase_max)
    def ralpha(self):
        return uniform(self.alpha_min,self.alpha_max)
    
    def rbeta(self):
        return uniform(self.beta_min,self.beta_max)
    
    def rrho(self):
        return uniform(self.rho_min,self.rho_max)
    
    def rpsi(self):
        return uniform(self.psi_min,self.psi_max)
    
    def rPeriod(self):
        return uniform(self.Pmin,self.Pmax)
    def rtheta(self):
        return uniform(self.theta_min,self.theta_max)

    def SenoidalPura(self):
        return self.A*sin(2*pi*self.f*self.t-self.rphi())
    
    def Sag(self):
        u = self.Interval()
        AFinal = self.A*(1-self.ralpha()*u)
        return AFinal*sin(2*pi*self.f*self.t - self.rphi()) 
    
    def Swell(self):
        u=self.Interval()
        AFinal=self.A*(1+self.rbeta()*u) #The points of the sag will be multiplied by the factor beta (their amplitude will increase). The rest, no.
        return AFinal*sin(2*pi*self.f*self.t-self.rphi()) #Final distorted sample
    
    def Interruption(self):
        u=self.Interval()
        AFinal=self.A*(1-self.rrho()*u) # The points of the sag will be multiplied by the factor rho (their amplitude will increase much). The rest, no.
        return AFinal*sin(2*pi*self.f*self.t-self.rphi()) #Final distorted sample

    def Impulse(self):
        # select transient interval (ocurrence point: tb-ta)
        taPeriod=self.rPeriod()
        ta=taPeriod*(1/self.f)
        pointts=round(taPeriod*(self.fs/self.f))

        u1=np.concatenate((np.zeros((pointts)), np.ones((self.PointsPerSignal-pointts))))# Step funcion, translated pointts to the right
        u2=np.concatenate((np.zeros(pointts+self.Onems), np.ones(self.PointsPerSignal-(pointts+self.Onems))))# Step funcion, translated (pointts+Onems) points to the right
        u=u1-u2 # Unitary function with all 0 except the interval in which the transient occurs (all 0 except tb-ta)
        # Final distorted sample
        exp_sub= (np.exp(-750*(self.t-ta)))-(np.exp(-344*(self.t-ta)))#;
        AFinal= self.A*self.rpsi()*exp_sub*u#;
        return -AFinal + self.A*sin(2*pi*self.f*self.t-self.rphi())#
    
    def Oscillatory(self):
        beta,fn,tau,theta=self.OscTranParam()#Select oscillatory transiente parameters
        u,t1=self.OscTranInterval()# Select oscillatory transient occurrence interval and starting and ending points
        return self.A* sin(2*pi*self.f*self.t-self.rphi())+self.A*beta*exp(-(self.t-t1)/tau)*sin(2*pi*fn*(self.t-t1)-theta)*u# Final distorted sample

    def Harmonics(self):
        alpha3,alpha5,theta1,theta3,theta5 = self.HarmParam()#Select harmonics parameters
        alpha7=uniform(self.alpha7_min,self.alpha7_max)#;
        theta7=uniform(self.theta_min,self.theta_max)#;
        return self.A*(self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5)+alpha7*sin(7*2*pi*self.f*self.t-theta7))#Final distorted sample


    def HarmonicSag(self):
        alpha3,alpha5,theta1,theta3,theta5=self.HarmParam()
        #  Random sag interval selection
        u = self.Interval()
        #  Final distorted sample
        AFinal=self.A*(1-self.ralpha()*u)
        return self.A*AFinal*(self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))

    def HarmonicSwell(self):
        alpha3,alpha5,theta1,theta3,theta5=self.HarmParam()
        #  Random swell interval selection
        u=self.Interval()
        #  Final distorted sample
        AFinal=self.A*(1+self.rbeta()*u)
        return self.A*AFinal*(self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))

    def Flicker(self):
        lbda,ff,phi=self.FlickerParam()
        #Final distorted sample
        return self.A*sin(2*pi*self.f*self.t-phi)*(1+lbda*sin(2*pi*ff*self.t))
    
    def FlickerSag(self):
        # Random amplitude, phase and frequency parameters selection
        lbda,ff,phi=self.FlickerParam()
        # Random sag interval selection
        u=self.Interval()
        # Final distorted sample
        return (self.A*sin(2*pi*self.f*self.t-phi))*(lbda*sin(2*pi*ff*self.t)+(1-self.ralpha()*u))


    def FlickerSwell(self):
        
        lbda,ff,phi = self.FlickerParam()
        # Random swell interval selection
        u = self.Interval()
        # Final distorted sample
        return (self.A*sin(2*pi*self.f*self.t-phi))*(lbda*sin(2*pi*ff*self.t)+(1+self.rbeta()*u))

    def OscillatoryTrasient(self):
        beta,fn,tau,theta=self.OscTranParam()
        # Random sag/osc transient interval selection
        u,utran,t1tran=self.OscTranInSagSwellInterval()
        # Final distorted sample
        return self.A*sin(2*pi*self.f*self.t-self.rphi())*(1-self.ralpha()*u)+self.A*beta*exp(-(self.t-t1tran)/tau)*sin(2*pi*fn*(self.t-t1tran)-theta)*utran

    def SagHarmonic(self):
        alpha3,alpha5,theta1,theta3,theta5=self.HarmParam()
        #  Random sag interval selection
        u=self.Interval()
        #  Final distorted sample
        return self.A*(self.alpha1*sin(2*pi*self.f*self.t-theta1)+(self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(-self.ralpha()*u))

    def SwellHarmonic(self):
        # Random selection of swell and harmonics parameters
        alpha3,alpha5,theta1,theta3,theta5=self.HarmParam()
        # Random swell interval selection
        u = self.Interval()
        # Final distorted sample
        return self.A*(self.alpha1*sin(2*pi*self.f*self.t-theta1)+(self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(self.rbeta()*u))

    def Notch(self):
        # Select number of notchs per period randomly
        notch_number=self.c[randint(0,len(self.c)-1)]
            # Select notch time occurrence
        td_max=(1/(notch_number*self.f))
        factor=1/(notch_number*self.f)
        ut=np.zeros(self.PointsPerSignal)
        k=uniform(self.k_min,self.k_max)
        tdminustc= uniform(self.tdminustc_min,self.tdminustc_max)
        td=uniform(self.td_min,td_max)
        tc=td-tdminustc

        while(tc<self.tc_min):
            td=uniform(self.td_min,td_max)
            tc=td-tdminustc

        for nn in range((self.n)*notch_number-1):

            pointt1=round((tc+factor*nn)*(self.fs)) # The period is converted into points
            pointt2=round((td+factor*nn)*(self.fs))
            u1=np.concatenate((np.zeros(pointt1) ,np.ones(self.PointsPerSignal-(pointt1))))
            u2=np.concatenate((np.zeros(pointt2),np.ones(self.PointsPerSignal-(pointt2))))
            u=k*(u1-u2) # Difference function, 0 all except the interval (tc+T*nn,td+T*nn)
            ut=ut+u # Final unitary function with all 0 except the notch occurrences
            
            # Final distorted sample
        return self.A*(sin(2*np.pi*self.f*self.t-self.rphi())-np.sign(np.sin(2*np.pi*self.f*self.t-self.rphi()))*ut)




    def HarmonicSagFlicker(self):
        # Selection of sag and harmonics parameters randomly
        alpha3,alpha5,theta1,theta3,theta5=self.HarmParam()
        # Selection of flicker parameters randomly
        lbda= uniform(self.lambda_min,self.lambda_max)
        ff = uniform(self.ff_min,self.ff_max)
        # Random sag interval selection
        u=self.Interval()
        # Final distorted sample
        AFinal = self.A*(1-self.ralpha()*u)
        return self.A*AFinal*(self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(1+lbda*sin(2*pi*ff*self.t))

    def HarmonicSwellFlicker(self):
        alpha3,alpha5,theta1,theta3,theta5=self.HarmParam()
        # Selection of flicker parameters randomly
        lbda= uniform(self.lambda_min,self.lambda_max)
        ff = uniform(self.ff_min,self.ff_max)
        #  Random swell interval selection
        u=self.Interval()
        #  Final distorted sample
        AFinal=self.A*(1+self.rbeta()*u)
        return self.A*AFinal*(self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(1+lbda*sin(2*pi*ff*self.t))

    def SagHarmonicFlicker(self):
        #Selection of sag and harmonics parameters randomly        
        alpha3,alpha5,theta1,theta3,theta5=self.HarmParam()
        # Selection of flicker parameters randomly
        lbda= uniform(self.lambda_min,self.lambda_max)
        ff = uniform(self.ff_min,self.ff_max)
        #  Random swell interval selection
        u=self.Interval()
        #  Final distorted sample

        return self.A*(self.alpha1*sin(2*pi*self.f*self.t-theta1)+(self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(-self.ralpha()*u)*(1+lbda*sin(2*pi*ff*self.t)))

    def SwellHarmonicFlicker(self):
        # Selection of swell and harmonics parameters randomly
        alpha3,alpha5,theta1,theta3,theta5=self.HarmParam()
        # Selection of flicker parameters randomly
        lbda= uniform(self.lambda_min,self.lambda_max)
        ff = uniform(self.ff_min,self.ff_max)
        #  Random swell interval selection
        u=self.Interval()
        #  Final distorted sample
        return self.A*(self.alpha1*sin(2*pi*self.f*self.t-theta1)+(self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(self.rbeta()*u)*(1+lbda*sin(2*pi*ff*self.t)))

    def SagHarmonicOscillatory(self):
        #Selection of sag and harmonics parameters randomly
        alpha3,alpha5,theta1,theta3,theta5=self.HarmParam()
        # Selection of oscillatory transient parametras randomly
        beta,fn,tau,theta=self.OscTranParam()
        # Random sag/osc transient interval selection
        u,utran,t1tran=self.OscTranInSagSwellInterval()
        # Final distorted sample
        return self.A*(self.alpha1*sin(2*pi*self.f*self.t-theta1)+(self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(-self.ralpha()*u)+beta*exp(-(self.t-t1tran)/tau)*sin(2*pi*fn*(self.t-t1tran)-theta)*utran)


    def SwellHarmonicOscillatory(self):
        # Selection of swell and harmonics parameters randomly
        alpha3,alpha5,theta1,theta3,theta5=self.HarmParam()
        #  Selection of oscillatory transient parametras randomly
        beta2,fn,tau,theta=self.OscTranParam()
        #  Random swell/osc transient interval selection
        u,utran,t1tran=self.OscTranInSagSwellInterval()
        #  Final distorted sample
        return self.A*(self.alpha1*sin(2*pi*self.f*self.t-theta1)+(self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(self.rbeta()*u)+beta2*exp(-(self.t-t1tran)/tau)*sin(2*pi*fn*(self.t-t1tran)-theta)*utran)


    def HarmonicSagOscillatory(self):

        # Selection of sag and harmonics parameters randomly
        alpha3,alpha5,theta1,theta3,theta5=self.HarmParam()
        # Random sag interval selection
        u=self.Interval()
        # Random oscillatory transient interval selection
        utran,t1tran  =self.OscTranInterval()
        # Selection of oscillatory transient parameters randomly
        beta,fn,tau,theta = self.OscTranParam()
        # Final distorted sample
        return self.A*((self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(1-self.ralpha()*u)+beta*exp(-(self.t-t1tran)/tau)*sin(2*pi*fn*(self.t-t1tran)-theta)*utran)

    def HarmonicSwellOscillatory(self):
        alpha3,alpha5,theta1,theta3,theta5 = self.HarmParam()
        # Random swell interval selection
        u=self.Interval()
        # Random oscillatory transient interval selection
        utran,t1tran =self.OscTranInterval()
        # Selection of oscillatory transient parameters randomly
        beta2,fn,tau,theta = self.OscTranParam()
        # Final distorted sample
        return self.A*((self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(1+self.rbeta()*u)+beta2*exp(-(self.t-t1tran)/tau)*sin(2*pi*fn*(self.t-t1tran)-theta)*utran)
    
    def HarmonicSagFlickerOscillatory(self):
        
        alpha3,alpha5,theta1,theta3,theta5=self.HarmParam()
        #  Random sag interval selection
        u = self.Interval()
        #  Flicker
        lbda= uniform(self.lambda_min,self.lambda_max)
        ff = uniform(self.ff_min,self.ff_max)

        utran,t1tran = self.OscTranInterval()
        #  Selection of oscillatory transient parameters randomly
        beta,fn,tau,theta = self.OscTranParam()
        #  Final distorted sample
        return self.A*(((self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(1-self.ralpha()*u))+beta*exp(-(self.t-t1tran)/tau)*sin(2*pi*fn*(self.t-t1tran)-theta)*utran)*(1+lbda*sin(2*pi*ff*self.t))



    def HarmonicSwellFlickerOscillatory(self):
        # Selection of swell and harmonics parameters randomly
        alpha3,alpha5,theta1,theta3,theta5=self.HarmParam()
        #  Random swell interval selection
        u=self.Interval()
        #  Flicker
        lbda= uniform(self.lambda_min,self.lambda_max)
        ff = uniform(self.ff_min,self.ff_max)

        utran,t1tran = self.OscTranInterval()
        #  Selection of oscillatory transient parameters randomly
        beta,fn,tau,theta = self.OscTranParam()
        #  Final distorted sample
        return self.A*(((self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(1+self.rbeta()*u))+beta*exp(-(self.t-t1tran)/tau)*sin(2*pi*fn*(self.t-t1tran)-theta)*utran)*(1+lbda*sin(2*pi*ff*self.t))

    def SagHarmonicFlickerOscillatory(self):

        alpha3,alpha5,theta1,theta3,theta5 = self.HarmParam()
        # Flicker
        lbda= uniform(self.lambda_min,self.lambda_max)
        ff = uniform(self.ff_min,self.ff_max)
        # Random sag/osc. transient interval selection
        u,utran,t1tran = self.OscTranInSagSwellInterval()
        # Selection of oscillatory transient parameters randomly
        beta,fn,tau,theta = self.OscTranParam()
        # Final distorted sample
        return self.A*(self.alpha1*sin(2*pi*self.f*self.t-theta1)+((self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(-self.ralpha()*u)+beta*exp(-(self.t-t1tran)/tau)*sin(2*pi*fn*(self.t-t1tran)-theta)*utran)*(1+lbda*sin(2*pi*ff*self.t)))


    def SwellHarmonicFlickerOscillatory(self):
    # Selection of swell and harmonics parameters randomly
        alpha3,alpha5,theta1,theta3,theta5 = self.HarmParam()
        #  Flicker
        lbda= uniform(self.lambda_min,self.lambda_max)
        ff = uniform(self.ff_min,self.ff_max)
        #  Random sag/osc. transient interval selection
        u,utran,t1tran = self.OscTranInSagSwellInterval()
        #  Selection of oscillatory transient parameters randomly
        beta2,fn,tau,theta = self.OscTranParam()
        #  Final distorted sample
        return self.A*(self.alpha1*sin(2*pi*self.f*self.t-theta1)+((self.alpha1*sin(2*pi*self.f*self.t-theta1)+alpha3*sin(3*2*pi*self.f*self.t-theta3)+alpha5*sin(5*2*pi*self.f*self.t-theta5))*(self.rbeta()*u)+beta2*exp(-(self.t-t1tran)/tau)*sin(2*pi*fn*(self.t-t1tran)-theta)*utran)*(1+lbda*sin(2*pi*ff*self.t)))

    def PQaleatorio(self,N):
        tamano=29
        Senales = np.zeros((N,tamano,self.PointsPerSignal))
        for v in range (N):
            #for k in range (tamano):
            Senales[v,0] = self.SenoidalPura()
            Senales[v,1] = self.Sag()
            Senales[v,2] = self.Swell()
            Senales[v,3] = self.Interruption()
            Senales[v,4] = self.Impulse()
            Senales[v,5] = self.Oscillatory()
            Senales[v,6] = self.Harmonics()
            Senales[v,7] = self.HarmonicSag()
            Senales[v,8] = self.HarmonicSwell()
            Senales[v,9] = self.Flicker()            
            Senales[v,10] = self.FlickerSag()
            Senales[v,11] = self.FlickerSwell()
            Senales[v,12] = self.SagHarmonicOscillatory()
            Senales[v,13] = self.SwellHarmonicOscillatory()
            Senales[v,14] = self.SagHarmonic()
            Senales[v,15] = self.SwellHarmonic()
            Senales[v,16] = self.Notch()
            Senales[v,17] = self.HarmonicSagFlicker()
            Senales[v,18] = self.HarmonicSwellFlicker()
            Senales[v,19] = self.SagHarmonicFlicker()
            Senales[v,20] = self.SwellHarmonicFlicker()            
            Senales[v,21] = self.SagHarmonicOscillatory()
            Senales[v,22] = self.SwellHarmonicOscillatory()
            Senales[v,23] = self.HarmonicSagOscillatory()
            Senales[v,24] = self.HarmonicSwellOscillatory()
            Senales[v,25] = self.HarmonicSagFlickerOscillatory()
            Senales[v,26] = self.HarmonicSwellFlickerOscillatory()
            Senales[v,27] = self.SagHarmonicFlickerOscillatory()
            Senales[v,28] = self.SwellHarmonicFlickerOscillatory() 

        return Senales
######################Auxiliary funcionts##########################################

    def OscTranInSagSwellInterval(self):
        # Random sag or swell interval selection
        if self.t1t2 == None:
            t1t2 = uniform(self.Pmin,self.Pmax) # Random number between periodMin-periodMax indicating how many periods the disturbace last
            pointst1t2 = round((t1t2)*(self.fs/self.f))   # The period is converted into points
        else:     
            pointst1t2 = round((self.t1t2 )*(self.fs/self.f))   # The period is converted into points
            border = round((self.PeriodLimit)*(self.fs/self.f))

        if self.Inicio == None:
            pointt1 = round((self.PointsPerSignal)*random()) # Random initial point of the disturbace selected randomly in the range 0-PointsPerSignal

            while (pointt1+pointst1t2)>self.PointsPerSignal-border:   # We check that the disturbace ends before the signal ends, otherwise a new starting point is generated
                pointt1 = round((self.PointsPerSignal)*random()) #Initial point of the disturbace selected randomly in the range 0-PointsPerSignal
        else:         
            pointt1=round((self.Inicio)*(self.fs/self.f))

        u1=np.concatenate((np.zeros(pointt1),np.ones(self.PointsPerSignal-pointt1)))
        u2=np.concatenate((np.zeros(pointt1+pointst1t2), np.ones(self.PointsPerSignal-(pointt1+pointst1t2))))
        u=u1-u2#
        # Random oscillatory transient interval selection
        pointst1t2tran= round(self.pointsfithpartperiod+(pointst1t2-self.pointsfithpartperiod)*random()) # Random number to select the duration of the oscillatory transient (it has to be whithin the sag duration)
        pointt1tran=round(pointt1+((pointt1+pointst1t2-pointst1t2tran)-pointt1)*random()) # Initial point of the oscillatory transient disturbace
        t1tran=pointt1tran*(1/self.fs) # Convert intial oscillatory transient point into time
        u1=np.concatenate((np.zeros(pointt1tran), np.ones(self.PointsPerSignal-pointt1tran))) # Step funcion, translated pointt1tran points to the right
        u2=np.concatenate((np.zeros(pointt1tran+pointst1t2tran), np.ones(self.PointsPerSignal-(pointt1tran+pointst1t2tran)))) # Step funcion, translated pointt1tran+pointst1t2tran points to the right
        
        utran=u2-u1 #  Difference function, 0 all except the interval (pointt1tran,pointt1tran+pointst1t2tran)

        return u,utran,t1tran

    def FlickerParam(self):
        lamb = uniform(self.lambda_min,self.lambda_max)
        ff = uniform(self.ff_min,self.ff_max)
        phi = uniform(self.phase_min,self.phase_min)
        return lamb,ff,phi

    def HarmParam(self):
        alpha3 = uniform(self.alpha3_min,self.alpha3_max)
        alpha5 = uniform(self.alpha5_min,self.alpha5_max)
        theta1 = uniform(self.theta_min,self.theta_max)
        theta3 = uniform(self.theta_min,self.theta_max)
        theta5 = uniform(self.theta_min,self.theta_max)
        return alpha3,alpha5,theta1,theta3,theta5
    
    def OscTranInterval(self):
        '''t1t2 = uniform(self.periodMinOT,self.periodMaxOT) #Random number between periodMinOT-periodMaxOT, indicating how many periods the disturbace last
        pointst1t2 = round(t1t2*(self.fs/self.f)) #The period is converted into points
        pointt1 = round((self.PointsPerSignal)*random()) #Initial point of the disturbace: random number between 0 and PointsPerSignal
        while (pointt1+pointst1t2)>self.PointsPerSignal:
            pointt1 = round((self.PointsPerSignal-0)*random())'''

        if self.t1t2 == None:
            t1t2 = uniform(self.periodMinOT,self.periodMaxOT) #Random number between periodMinOT-periodMaxOT, indicating how many periods the disturbace last
            pointst1t2 = round((t1t2)*(self.fs/self.f))   # The period is converted into points
        else:     
            pointst1t2 = round((self.t1t2 )*(self.fs/self.f))   # The period is converted into points
            border = round((self.PeriodLimit)*(self.fs/self.f))

        if self.Inicio == None:
            pointt1 = round((self.PointsPerSignal)*random()) # Random initial point of the disturbace selected randomly in the range 0-PointsPerSignal

            while (pointt1+pointst1t2)>self.PointsPerSignal-border:   # We check that the disturbace ends before the signal ends, otherwise a new starting point is generated
                pointt1 = round((self.PointsPerSignal)*random()) #Initial point of the disturbace selected randomly in the range 0-PointsPerSignal
        else:         
            pointt1=round((self.Inicio)*(self.fs/self.f))    

        t1 = pointt1*(1/self.fs) #convert pointt1 into time
        u1 = np.concatenate([np.zeros((pointt1)), np.ones((self.PointsPerSignal-pointt1))]); #Step funcion, translated pointt1 points to the right
        u2 = np.concatenate([np.zeros((pointt1+pointst1t2)),np.ones((self.PointsPerSignal-(pointt1+pointst1t2)))]) #Step funcion, translated pointt1+pointst1t2 points to the right
        u = u1-u2 # Difference function, 0 all except the interval (pointt1,pointt1+pointst1t2)

        return u,t1

    def OscTranParam(self):
        beta = self.rbeta()
        fn= uniform(self.fn_min,self.fn_max)
        tau= uniform(self.tau_min,self.tau_max)
        theta= self.rtheta()
        return beta,fn,tau,theta    
    
    
    def Interval(self):

        if self.t1t2 == None:
            t1t2 = uniform(self.Pmin,self.Pmax) # Random number between periodMin-periodMax indicating how many periods the disturbace last
            pointst1t2 = round((t1t2)*(self.fs/self.f))   # The period is converted into points
        else:     
            pointst1t2 = round((self.t1t2 )*(self.fs/self.f))   # The period is converted into points
            border = round((self.PeriodLimit)*(self.fs/self.f))

        if self.Inicio == None:
            pointt1 = round((self.PointsPerSignal)*random()) # Random initial point of the disturbace selected randomly in the range 0-PointsPerSignal

            while (pointt1+pointst1t2)>self.PointsPerSignal-border:   # We check that the disturbace ends before the signal ends, otherwise a new starting point is generated
                pointt1 = round((self.PointsPerSignal)*random()) #Initial point of the disturbace selected randomly in the range 0-PointsPerSignal
        else:         
            pointt1=round((self.Inicio)*(self.fs/self.f))
            
        u1 = np.concatenate((np.zeros(pointt1), np.ones(self.PointsPerSignal-pointt1))) #Step funcion, translated pointt1 points to the right
        u2 = np.concatenate((np.zeros(pointt1+pointst1t2),np.ones(self.PointsPerSignal-(pointt1+pointst1t2)))) #Step funcion, translated pointt1+pointst1t2 points to the right
        u = u1-u2 # Difference function, 0 all except the interval (pointt1,pointt1+pointst1t2)

        return u
    def FramesAnimacion(self):



        return self.SenoidalPura()






#prueba=PQ(Cicles=17,PeriodoDisturbio=5,FinalDisturbio=2,InicioDisturbio=5)

#plt.plot(prueba.t,prueba.PQaleatorio(300)[1,28])
#plt.show()

