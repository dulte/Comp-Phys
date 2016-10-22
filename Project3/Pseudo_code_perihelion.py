previous_r =0
previous_previous_r=0

for k in range(0, time, N):
    r_vec = mercury_position - sun_position
    r = r_vec.length();
    if previous_r < r and previous_r < previous_previous_r:
        theta = arctan2(r_vec.y,r_vec.x);  //Save as perihelion passing
