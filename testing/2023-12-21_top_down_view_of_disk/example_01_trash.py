import matplotlib.pyplot as plt
import numpy as np

# Function to plot a circle representing a particle in the protoplanetary disk
def plot_particle(ax, radius, center, color='b'):
    circle = plt.Circle(center, radius, color=color, fill=False)
    ax.add_patch(circle)

# Function to create a top-down view of a protoplanetary disk
def plot_protoplanetary_disk(num_particles=100, disk_radius=1.0):
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot central star
    plot_particle(ax, 0.1, (0, 0), color='yellow')  # Assuming the star is yellow

    # Plot protoplanetary disk particles
    for _ in range(num_particles):
        radius = np.random.uniform(0.01, 0.1)  # Random particle size
        distance = np.random.uniform(0, disk_radius)
        angle = np.random.uniform(0, 2 * np.pi)
        x = distance * np.cos(angle)
        y = distance * np.sin(angle)
        plot_particle(ax, radius, (x, y))

    # Set aspect ratio to be equal
    ax.set_aspect('equal', adjustable='box')

    # Set plot limits
    ax.set_xlim(-disk_radius, disk_radius)
    ax.set_ylim(-disk_radius, disk_radius)

    # Set labels and title
    ax.set_xlabel('X-axis (AU)')
    ax.set_ylabel('Y-axis (AU)')
    ax.set_title('Protoplanetary Disk - Top-Down View')

    # Show the plot
    plt.show()

# Example usage
plot_protoplanetary_disk(num_particles=100, disk_radius=1.0)
