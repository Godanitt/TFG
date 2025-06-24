import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Ajuste global del grosor de las líneas
mpl.rcParams['lines.linewidth'] = 2.5  # o el valor que desees


Energias=np.array([0,6.2,9.8,16.2,19.8,21.4,26.4,31.8,34.4,36.6,39.2,49.8,52.2,55.2,57.6,59.6,72.2,75.0,77.8,79.6,81.6])
Degeneracion=np.array([2,4,2,6,2,4,8,4,6,2,10,8,6,4,2,12,10,8,6,4,2])
Nombre=["1s1/2  [0MeV]","1p3/2 [6.2MeV]","1p1/2 [9.8MeV]","1d5/2 [16.2MeV]",    "2s1/2 [19.8MeV]",  "1d3/2 [21.4MeV]",    "1f7/2 [26.4MeV]", 
        "2p3/2 [31.8MeV]",   "1f5/2 [34.4MeV]",    "2p1/2 [36.6MeV]",   "1g9/2 [39.2MeV]",    "1g7/2 [49.8 MeV]", 
        "2d5/2 [52.2MeV]",    "2d3/2 [55.2MeV]",    "3s1/2 [57.6MeV]",    "1h11/2 [59.6MeV]",    "1h9/2 [72.2MeV]",  
        "2f7/2 [75.0MeV]",    "2f5/2[77.8MeV]",     "3p3/2 [79.6MeV]",    "3p1/2 [81.6 MeV]"]

Nombre_sin_energia = [
    "$1s_{1/2}$", "$1p_{3/2}$", "$1p_{1/2}$", "$1d_{5/2}$", "$2s_{1/2}$", "$1d_{3/2}$", "1f7/2",
    "2p3/2", "1f5/2", "2p1/2", "1g9/2", "1g7/2", "2d5/2", "2d3/2",
    "3s1/2", "1h11/2", "1h9/2", "2f7/2", "2f5/2", "3p3/2", "3p1/2"
]

Nombre_sin_energia2 = [
    "$1s_{1/2}$", "$1p_{3/2}$","$2s_{1/2}$", "$1p_{1/2}$", "$1d_{5/2}$",  "$1d_{3/2}$", "1f7/2",
    "2p3/2", "1f5/2", "2p1/2", "1g9/2", "1g7/2", "2d5/2", "2d3/2",
    "3s1/2", "1h11/2", "1h9/2", "2f7/2", "2f5/2", "3p3/2", "3p1/2"
]






def Modelo_Capas(A,Z,nombre):
    N=A-Z   
    fig,ax= plt.subplots()
    
    for i in range(len(Degeneracion)):
        
        if sum(Degeneracion[0:i+1])>N:
                
                
            plt.plot([1,5],[Energias[i],Energias[i]],color="black")
            
            
            plt.plot([1,5],[Energias[i+1],Energias[i+1]],color="black",linestyle="--")
            plt.plot([-5,-1],[Energias[i+1],Energias[i+1]],color="black",linestyle="--")
            plt.annotate(nombre[i+1],xy=(-0.35,Energias[i+1]-0.2))
            
            plt.plot([1,5],[Energias[i+2],Energias[i+2]],color="black",linestyle="--")
            plt.plot([-5,-1],[Energias[i+2],Energias[i+2]],color="black",linestyle="--")
            plt.annotate(nombre[i+2],xy=(-0.35,Energias[i+2]-0.2))
            
            aux=N-sum(Degeneracion[0:i])
            nucleon=np.linspace(1,5,2+aux)
            nucleon=nucleon[1:2+aux-1]
            n=len(nucleon)
            plt.scatter(nucleon,[Energias[i]]*n,color="grey",
                s=200,                # Tamaño del marker (área en puntos^2)
                edgecolors='black',   # Color del borde
                linewidths=1.5,         # Grosor del borde
                marker='o')           # Tipo de marker
            
            break
        
        plt.plot([1,5],[Energias[i],Energias[i]],color="black")
        plt.plot([-5,-1],[Energias[i],Energias[i]],color="black",linestyle="--")
        plt.annotate(nombre[i],xy=(-0.35,Energias[i]-0.2))
        
        nucleon=np.linspace(1,5,2+Degeneracion[i])
        nucleon=nucleon[1:2+Degeneracion[i]-1]
        n=len(nucleon)
        plt.scatter(nucleon,[Energias[i]]*n,color="grey",
            s=200,                # Tamaño del marker (área en puntos^2)
            edgecolors='black',   # Color del borde
            linewidths=1.5,         # Grosor del borde
            marker='o')           # Tipo de marker
        
    
           
    for i in range(len(Degeneracion)):
        if sum(Degeneracion[0:i+1])>Z:
                
            plt.plot([-5,-1],[Energias[i],Energias[i]],color="black")
            
            plt.plot([1,5],[Energias[i+1],Energias[i+1]],color="black",linestyle="--")
            plt.plot([-5,-1],[Energias[i+1],Energias[i+1]],color="black",linestyle="--")
            plt.annotate(nombre[i+1],xy=(-0.35,Energias[i+1]-0.2))
            
            plt.plot([1,5],[Energias[i+2],Energias[i+2]],color="black",linestyle="--")
            plt.plot([-5,-1],[Energias[i+2],Energias[i+2]],color="black",linestyle="--")
            plt.annotate(nombre[i+2],xy=(-0.35,Energias[i+2]-0.2))
            
            aux=Z-sum(Degeneracion[0:i])
            nucleon=np.linspace(-5,-1,2+aux)
            nucleon=nucleon[1:2+aux-1]
            n=len(nucleon)
            plt.scatter(nucleon,[Energias[i]]*n,color="red",
                s=200,                # Tamaño del marker (área en puntos^2)
                edgecolors='black',   # Color del borde
                linewidths=1.5,         # Grosor del borde
                marker='o')           # Tipo de marker
            
    
            break
        
        plt.plot([-5,-1],[Energias[i],Energias[i]],color="black")
        plt.plot([1,5],[Energias[i],Energias[i]],color="black",linestyle="--")
        
        nucleon=np.linspace(-5,-1,2+Degeneracion[i])
        nucleon=nucleon[1:2+Degeneracion[i]-1]
        n=len(nucleon)
        plt.scatter(nucleon,[Energias[i]]*n,color="red",
            s=200,                # Tamaño del marker (área en puntos^2)
            edgecolors='black',   # Color del borde
            linewidths=1.5,         # Grosor del borde
            marker='o')           # Tipo de marker
        
    for i in range(len(Degeneracion)):
        
        if sum(Degeneracion[0:i+1])>N:
                
            
            aux=N-sum(Degeneracion[0:i])
            nucleon=np.linspace(1,5,2+aux)
            nucleon=nucleon[1:2+aux-1]
            n=len(nucleon)
            plt.scatter(nucleon,[Energias[i]]*n,color="grey",
                s=200,                # Tamaño del marker (área en puntos^2)
                edgecolors='black',   # Color del borde
                linewidths=1.5,         # Grosor del borde
                marker='o',zorder=5)           # Tipo de marker
            
            break
        
        
        nucleon=np.linspace(1,5,2+Degeneracion[i])
        nucleon=nucleon[1:2+Degeneracion[i]-1]
        n=len(nucleon)
        plt.scatter(nucleon,[Energias[i]]*n,color="grey",
            s=200,                # Tamaño del marker (área en puntos^2)
            edgecolors='black',   # Color del borde
            linewidths=1.5,         # Grosor del borde
            marker='o',zorder=5)           # Tipo de marker
        
    
           
    for i in range(len(Degeneracion)):
        if sum(Degeneracion[0:i+1])>Z:
                
            aux=Z-sum(Degeneracion[0:i])
            nucleon=np.linspace(-5,-1,2+aux)
            nucleon=nucleon[1:2+aux-1]
            n=len(nucleon)
            plt.scatter(nucleon,[Energias[i]]*n,color="red",
                s=200,                # Tamaño del marker (área en puntos^2)
                edgecolors='black',   # Color del borde
                linewidths=1.5,         # Grosor del borde
                marker='o',zorder=5)           # Tipo de marker
            
    
            break
        
        nucleon=np.linspace(-5,-1,2+Degeneracion[i])
        nucleon=nucleon[1:2+Degeneracion[i]-1]
        n=len(nucleon)
        plt.scatter(nucleon,[Energias[i]]*n,color="red",
            s=300,                # Tamaño del marker (área en puntos^2)
            edgecolors='black',   # Color del borde
            linewidths=1.5,         # Grosor del borde
            marker='o',zorder=5)           # Tipo de marker
    
    
    # Eliminar todos los bordes (caja)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    
    plt.annotate("Protones",xy=(-3.5,-3))
    plt.annotate("Neutrones",xy=(2.2,-3))
    ax.set_ylim(-4,21)
    ax.set_xlim(-5,5)
    
    
    return fig,ax

fig,ax=Modelo_Capas(10,3,Nombre_sin_energia)
fig.savefig("/home/daniel/GitHub/TFG/Memoria/Imagenes/Capas_10Li.pdf",bbox_inches="tight")

fig,ax=Modelo_Capas(10,3,Nombre_sin_energia2)
fig.savefig("/home/daniel/GitHub/TFG/Memoria/Imagenes/Capas_10Li2.pdf",bbox_inches="tight")