import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Ajuste global del grosor de las líneas
mpl.rcParams['lines.linewidth'] = 2.5  # o el valor que desees


Energias1=np.array([0,6.2,9.8,16.2,19.8,21.4,26.4,31.8,34.4,36.6,39.2,49.8,52.2,55.2,57.6,59.6,72.2,75.0,77.8,79.6,81.6])
Energias2=np.array([0,6.2,9.2,12.5,17.0,21.4,26.4,31.8,34.4,36.6,39.2,49.8,52.2,55.2,57.6,59.6,72.2,75.0,77.8,79.6,81.6])


EnergiasMagicos=np.array([3.1,13.0])
NombreMagicos=np.array(["2","8"])

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




# Configuración global del tamaño de fuente
mpl.rcParams['font.size'] = 20  # Cambia este valor a lo que necesites


def Modelo_Capas(A,Z,nombre,Energias=Energias1,flag=False):
    N=A-Z   
    fig,ax= plt.subplots()
    
    for i in range(len(Degeneracion)):
        
        if sum(Degeneracion[0:i+1])>N:
                
                
            plt.plot([1,5],[Energias[i],Energias[i]],color="black")
            plt.annotate(nombre[i],xy=(-0.75,Energias[i]-0.2))

            
            plt.plot([1,5],[Energias[i+1],Energias[i+1]],color="black",linestyle="--")
            plt.plot([-5,-1],[Energias[i+1],Energias[i+1]],color="black",linestyle="--")
            plt.annotate(nombre[i+1],xy=(-0.75,Energias[i+1]-0.2))
            
            plt.plot([1,5],[Energias[i+2],Energias[i+2]],color="black",linestyle="--")
            plt.plot([-5,-1],[Energias[i+2],Energias[i+2]],color="black",linestyle="--")
            plt.annotate(nombre[i+2],xy=(-0.75,Energias[i+2]-0.2))
            
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
        plt.annotate(nombre[i],xy=(-0.75,Energias[i]-0.2))
        
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
            
            plt.plot([1,5],[Energias[i+2],Energias[i+2]],color="black",linestyle="--")
            plt.plot([-5,-1],[Energias[i+2],Energias[i+2]],color="black",linestyle="--")
            
            aux=Z-sum(Degeneracion[0:i])
            nucleon=np.linspace(-5,-1,2+aux)
            nucleon=nucleon[1:2+aux-1]
            n=len(nucleon)
            plt.scatter(nucleon,[Energias[i]]*n,color="red",
                s=300,                # Tamaño del marker (área en puntos^2)
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
            s=300,                # Tamaño del marker (área en puntos^2)
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
                s=300,                # Tamaño del marker (área en puntos^2)
                edgecolors='black',   # Color del borde
                linewidths=1.5,         # Grosor del borde
                marker='o',zorder=5)           # Tipo de marker
            
            break
        
        
        nucleon=np.linspace(1,5,2+Degeneracion[i])
        nucleon=nucleon[1:2+Degeneracion[i]-1]
        n=len(nucleon)
        plt.scatter(nucleon,[Energias[i]]*n,color="grey",
            s=300,                # Tamaño del marker (área en puntos^2)
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
                s=300,                # Tamaño del marker (área en puntos^2)
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
    
    if (flag==True):
        ax.annotate(
                '',
                xy=(2,9.2), xycoords='data',
                xytext=(2,19.8), textcoords='data',
                arrowprops=dict(
                    arrowstyle='->,head_width=0.8,head_length=1',
                    connectionstyle='arc3,rad=0.3',
                    linewidth=5,
                    color=(0.2, 0.2, 0.2, 0.5)  # gris oscuro y 50% opacidad
                ),
                fontsize=12,
                color='black'
        )
        
    plt.annotate("Protones",xy=(-4.2,-3))
    plt.annotate("Neutrones",xy=(1.4,-3))
    ax.set_ylim(-4,21)
    ax.set_xlim(-5,5)
    
    
    return fig,ax

fig,ax=Modelo_Capas(10,3,Nombre_sin_energia)
for i in range(len(NombreMagicos)):
        ax.annotate(
                NombreMagicos[i],
                xy=(-3.0,EnergiasMagicos[i]), 
                fontsize=16,
                color='black',     ha="center", va="center",
                bbox=dict(boxstyle="circle",
                      fc="lightblue", ec="steelblue", lw=2)
        )
        ax.annotate(
                NombreMagicos[i],
                xy=(3,EnergiasMagicos[i]), 
                fontsize=16,
                color='black',   
                ha="center", va="center",
                bbox=dict(boxstyle="circle",
                      fc="lightblue", ec="steelblue", lw=2)  
        )
        
fig.savefig("/home/daniel/GitHub/TFG/Memoria/Imagenes/Capas_10Li.pdf",bbox_inches="tight")

fig,ax=Modelo_Capas(10,3,Nombre_sin_energia2,Energias=Energias2,flag=True)
i=5
ax.plot([1,5],[19.8,19.8],color="black",linestyle="--",alpha=0.1)
ax.plot([-5,-1],[19.8,19.8],color="black",linestyle="--",alpha=0.1)
ax.annotate("$2s_{1/2}$",xy=(-0.75,19.8-0.2),alpha=0.3)
fig.savefig("/home/daniel/GitHub/TFG/Memoria/Imagenes/Capas_10Li2.pdf",bbox_inches="tight")