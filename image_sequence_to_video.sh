rm *.webm
ffmpeg -i log10_density_%d.png -codec:v libvpx "wd_nova_log10_density_`git rev-parse HEAD`.webm"
ffmpeg -i log10_pressure_%d.png -codec:v libvpx "wd_nova_log10_pressure_`git rev-parse HEAD`.webm"
ffmpeg -i log10_temperature_%d.png -codec:v libvpx "wd_nova_log10_temperature_`git rev-parse HEAD`.webm"
ffmpeg -i x_velocity_%d.png -codec:v libvpx "wd_nova_x_velocity_`git rev-parse HEAD`.webm"
ffmpeg -i y_velocity_%d.png -codec:v libvpx "wd_nova_y_velocity_`git rev-parse HEAD`.webm"
