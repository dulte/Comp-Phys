infile=open('Planet_initial_data.txt')
for line in infile: #Read in planet data
	initial_position, initial_velocity, mass=line.split()
	if name != sun: #To ensure that COM is at the origin with 0 velocity
		system.create_particle(initial_position, initial_velocity, mass_of_planet)

total_position=sum(particle_mass*particle_position)
total_momentum=sum(particle_mass*particle_velocity)
sun_position=-total_position
sun_velocity=-total_momentum/mass_of_the_Sun

system.create_particle(sun_position, sun_velocity, mass_of_the_sun) #Create the sun

