
import matplotlib.animation as animation
import numpy as np
import matplotlib.pyplot as plt

from PQmodel import PQ

prueba=PQ(Cicles=500,PeriodoDisturbio=20,FinalDisturbio=2,InicioDisturbio=100)


# y=prueba.SenoidalPura()
# tipodesenal='Senoidal Pura'
# y=prueba.Sag()
# tipodesenal='Sag'
# y=prueba.Swell()
# tipodesenal='Swell'
# y=prueba.Interruption()
# tipodesenal='Interruption'
# y=prueba.Impulse()
# tipodesenal='Impulse'
# y=prueba.Oscillatory()
# tipodesenal='Oscillatory'
# y=prueba.Harmonics()
# tipodesenal='Harmonics'
# y=prueba.HarmonicSag()
# tipodesenal='HarmonicSag'
# y=prueba.HarmonicSwell()
# tipodesenal='HarmonicSwell'
# y=prueba.Flicker()     
# tipodesenal='Flicker'       
# y=prueba.FlickerSag()
# tipodesenal='FlickerSag'
# y=prueba.FlickerSwell()
# tipodesenal='FlickerSwell'
# y=prueba.SagHarmonicOscillatory()
# tipodesenal='SagHarmonicOscillatory'
# y=prueba.SwellHarmonicOscillatory()
# tipodesenal='SwellHarmonicOscillatory'
# y=prueba.SagHarmonic()
# tipodesenal='SagHarmonic'
# y=prueba.SwellHarmonic()
# tipodesenal='SwellHarmonic'
# y=prueba.Notch()
# tipodesenal='Notch'
# y=prueba.HarmonicSagFlicker()
# tipodesenal='HarmonicSagFlicker'
# y=prueba.HarmonicSwellFlicker()
# tipodesenal='HarmonicSwellFlicker'
# y=prueba.SagHarmonicFlicker()
# tipodesenal='SagHarmonicFlicker'
# y=prueba.SwellHarmonicFlicker()
# tipodesenal='SwellHarmonicFlicker'            
# y=prueba.SagHarmonicOscillatory()
# tipodesenal='SagHarmonicOscillatory'
# y=prueba.SwellHarmonicOscillatory()
# tipodesenal='SwellHarmonicOscillatory'
# y=prueba.HarmonicSagOscillatory()
# tipodesenal='HarmonicSagOscillatory'
# y=prueba.HarmonicSwellOscillatory()
# tipodesenal='HarmonicSwellOscillatory'
# y=prueba.HarmonicSagFlickerOscillatory()
# tipodesenal='HarmonicSagFlickerOscillatory'
# y=prueba.HarmonicSwellFlickerOscillatory()
# tipodesenal='HarmonicSwellFlickerOscillatory'
# y=prueba.SagHarmonicFlickerOscillatory()
# tipodesenal='SagHarmonicFlickerOscillatory'
y=prueba.SwellHarmonicFlickerOscillatory()
tipodesenal='SwellHarmonicFlickerOscillatory'

fig = plt.figure()

ax=fig.gca()
k=158*7


font = {
    'weight': 'normal',
    'size'  :  20,
    'color': 'lightgray'
}



def actualizar (i):
    ax.clear()
    if i*k < len(y) and (i-5) >= 0:
        #ax.set_ylim(0,2)
        ax.set_ylim(-2,2)    
        ax.plot(y[(i-5)*k:i*k])
        ax.text(0.95, 0.2, tipodesenal,
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes,
            fontdict=font)

    plt.title(tipodesenal)

anim=animation.FuncAnimation(fig=fig, func=actualizar,frames=100,interval = 270)
plt.show()