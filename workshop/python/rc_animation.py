import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

import seaborn as sns

###############################################

def CalculatePosition(radius,velocity,time,dt):
    """Calculates position of an object around a circle after some time with interval dt."""
    
    # Initial conditions
    theta = 0
    xini = radius * np.cos(theta)
    yini = radius * np.sin(theta)
    t = 0
    
    # Store positions and time
    xposition = [xini]
    yposition = [yini]
    storedtime = [t]
    
    # Calculate positions
    while t < time:
        t += dt
        x = radius * np.cos((velocity/radius)*t)
        y = radius * np.sin((velocity/radius)*t)
        
        xposition.append(x)
        yposition.append(y)
        storedtime.append(t)
    
    # Make an array
    xposition = np.array(xposition)
    yposition = np.array(yposition)
    storedtime = np.array(storedtime)
    
    return xposition,yposition,storedtime

###############################################

def MultiplePositions(radius,velocity,time,dt):
    """Calculate positions of multiple objects around a circle."""
    
    # Stop the calculation when the outermost point takes a whole revolution
    # Outermost point position
    outerposition = radius[len(radius)-1]
    
    # Calculate the positions of outermost point:
    xouter = CalculatePosition(radius[len(radius)-1],velocity[len(radius)-1],time,dt)[0]
    
    # Circumference of the outer circle
    circouter = 2*np.pi*outerposition
    
    # Distance the outer object traveled
    distance = 0
    istop = 0
        
    # New time
    storedtime = CalculatePosition(radius[0],velocity[0],time,dt)[2]
    for t in storedtime:
        if distance < circouter:
            distance = velocity[len(radius)-1]*t
        else:
            istop = np.where(storedtime == t)[0]     # find index of the numpy array 
            break
                   
    istop = int(istop)
    newstoredtime = storedtime[:istop-1]
            
    xmultiple = []
    ymultiple = []

    for i in range(len(radius)):
        x = CalculatePosition(radius[i],velocity[i],time,dt)[0]
        y = CalculatePosition(radius[i],velocity[i],time,dt)[1]
        xmultiple.append(x)
        ymultiple.append(y)

    xmultiple = np.array(xmultiple)
    ymultiple = np.array(ymultiple)
    
    return xmultiple, ymultiple, newstoredtime

###############################################

def PlotRotationCurve(radius,velocity,title,
                      xlabel='Radius (km)',ylabel='Velocity (km/s)',
                      xlim=1,
                      ylim=0.1):
    
    # Convert title to string
    title = str(title)
    
    # Plot
    fig = plt.figure(figsize=(6,6))
    ax = plt.axes()
    
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    ax.set_xlabel(xlabel,color='white')
    ax.set_ylabel(ylabel,color='white')
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    
    plt.title(title,color='white',fontsize='14')
    plt.scatter(radius,velocity,color='khaki')
    plt.plot(radius,velocity,color='white')
    plt.xlim(0,np.max(radius)+xlim)
    plt.ylim(0,np.max(velocity)+ylim)
    plt.show()
    
###############################################

def MakeAnimation(radius,velocity,time,dt,filename,title,
                 xlim=1,ylim=1,
                 size=False,masses=None):
    
    # Extract x and y positions, and time
    xpositions = MultiplePositions(radius,velocity,time,dt)[0]
    ypositions = MultiplePositions(radius,velocity,time,dt)[1]
    #storedtimes = CalculatePosition(radius,velocity,time,dt)[2]
    storedtimes = MultiplePositions(radius,velocity,time,dt)[2]
    
    # Sizes of dots based on masses
    if size == True:   # use an array of masses as an input for sizes
        area = [s * 5e-25 for s in masses]
    if size == False:  # use a default size
        area = 100
                      
    # Create a movie write object, set frame rate
    writer = ani.FFMpegWriter(fps=25)

    # Create a figure, 8"x8" in size
    fig = plt.figure(figsize=(8,8))
    
    # Change background color of the plot
    fig.patch.set_facecolor('black')
    ax = plt.axes()
    ax.set_facecolor('black')
    ax.set_facecolor('black')
    ax.set_xlabel('x (km)',color='white')
    ax.set_ylabel('y (km)',color='white')
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')
  
    # Convert filename and title to string
    filename = str(filename)
    title = str(title)

    # Set things up to save frames to a movie:
    #   fig = the figure the writer will record from
    with writer.saving(fig, filename, 100):

    # Loop
        i = 0            # start counter
        for t in storedtimes:
            plt.cla()
            for r in radius:
                    circle = plt.Circle((0, 0), r, color='white', fill=False)
                    plt.gca().add_patch(circle)
            plt.scatter(xpositions[:,i],ypositions[:,i],s=area,color='khaki')
            plt.suptitle(title,color='white',fontsize='18')
            plt.title('{:.1e} seconds'.format(t),color='white')
            plt.xlim(-np.max(radius)-xlim,np.max(radius)+xlim)
            plt.ylim(-np.max(radius)-1,np.max(radius)+ylim)
            plt.xlabel('x (km)',color='white',fontsize='14')
            plt.ylabel('y (km)',color='white',fontsize='14')
            ax.tick_params(axis='x', colors='white')
            ax.tick_params(axis='y', colors='white')
            i += 1
            
            # Save the current plot as a movie frame
            writer.grab_frame()
            
###############################################

# create a figure and axes
fig = plt.figure(figsize=(12,5))
ax1 = plt.subplot(1,2,1)   
ax2 = plt.subplot(1,2,2)

# set up the subplots as needed
ax1.set_xlim(( 0, 2))            
ax1.set_ylim((-2, 2))
ax1.set_xlabel('Time')
ax1.set_ylabel('Magnitude')

ax2.set_xlim((-2,2))
ax2.set_ylim((-2,2))
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_title('Phase Plane')

# create objects that will change in the animation. These are
# initially empty, and will be given new values for each frame
# in the animation.
txt_title = ax1.set_title('')
line1, = ax1.plot([], [], 'b', lw=2)     # ax.plot returns a list of 2D line objects
line2, = ax1.plot([], [], 'r', lw=2)
pt1, = ax2.plot([], [], 'g.', ms=20)
line3, = ax2.plot([], [], 'y', lw=2)

ax1.legend(['sin','cos']);

# animation function. This is called sequentially
def drawframe(n):
    x = np.linspace(0, 2, 1000)
    y1 = np.sin(2 * np.pi * (x - 0.01 * n))
    y2 = np.cos(2 * np.pi * (x - 0.01 * n))
    line1.set_data(x, y1)
    line2.set_data(x, y2)
    line3.set_data(y1[0:50],y2[0:50])
    pt1.set_data(y1[0],y2[0])
    txt_title.set_text('Frame = {0:4d}'.format(n))
    return (line1,line2)



from matplotlib import animation

# blit=True re-draws only the parts that have changed.
anim = animation.FuncAnimation(fig, drawframe, frames=100, interval=20, blit=True)



