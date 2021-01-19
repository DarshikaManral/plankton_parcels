import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4_3D
from datetime import timedelta, datetime

# zooplankton_drift kernel: 
# defines the vertical movement of a basic plankton particle that performs Diel Vertical Migration
def zooplankton_drift(particle, fieldset, time):    
    
    migration_speed = 0.08     #(m/s)
 
    min_depth = 45      #(m)
    max_depth = 400     #(m) Check this

    # maximum displacement in one dt=0.08 * 300= 24 m (in 5 minutes)
    # https://aslopubs.onlinelibrary.wiley.com/doi/pdf/10.1002/lno.10219   (FIGURE 6)
   
    max_displacement= migration_speed * particle.dt
    
    #TODO
    # day and night timings are hard coded- works well for equatorial waters
    # maybe use a function on latitude and time to calculate the daylight and night timings 
    # error handling
    
    # day from 6 am to 6 pm 
    # night from 6 pm to 6 am 

    # Extract the current hour of the day from the origin timestamp from velocity files(data available from)
    # and time since the origin time.
    
    total_seconds=fieldset.start_time + time
    total_hour=total_seconds/3600
    current_hour= math.fmod(total_hour, 24)
        
    if current_hour>=6 and current_hour < 18: #day
        if((particle.depth + max_displacement)>max_depth):
            particle.depth=max_depth
        else:
            particle.depth+=max_displacement
            
    else:   #night
        if (particle.depth - max_displacement)<min_depth:
            particle.depth=min_depth
        else:
            particle.depth-=max_displacement
    
    
data_path='/data/oceanparcels/input_data/NEMO-MEDUSA/ORCA025-N006/'

mesh_mask = data_path + 'domain/coordinates.nc'

ufiles = sorted(glob(data_path + 'means/ORCA025-N06_201501*d05U.nc'))
vfiles = sorted(glob(data_path + 'means/ORCA025-N06_201501*d05V.nc'))
wfiles = sorted(glob(data_path + 'means/ORCA025-N06_201501*d05W.nc'))

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles},
             'W': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': wfiles}
            }

variables = {'U': 'uo',
             'V': 'vo',
             'W': 'wo'}

dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
              'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
              'W': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}
             }

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, chunksize='auto')


# extract the origin timestamp from one of the velocity files
# first timestamp from the first file to be read by the model 
u_temp = nc.Dataset(ufiles[0])
ticks = u_temp['time_counter'][:][0]
time_zero = datetime(1900, 1, 1) + timedelta(seconds=ticks)
print(time_zero)
time_zero_totalseconds = time_zero.hour * 60 * 60 + time_zero.minute * 60 + time_zero.second
print("start time:  ", time_zero_totalseconds)

# to make this value available to the custom kernel during execution add it to fieldset
fieldset.add_constant('start_time',time_zero_totalseconds)

pset = ParticleSet.from_line(fieldset=fieldset, 
                             size=10, 
                             pclass=JITParticle,
                             start=(-47.07351,  1.50464), 
                             finish=(-42.23952, 1.50464), 
                             time=datetime(2015, 1, 3, 12, 0, 0), 
#                              repeatdt= timedelta(days=1),
                             depth=400) # in m; since start time is noon, if it starts at night, depth = 45 m

# two day simulation therefore, outputdt is small, else per day.

output_file_path="/scratch/manra003/Plankton_withDVM_short_run_2D.nc"
output_file = pset.ParticleFile(name=output_file_path, outputdt=300)

kernels = pset.Kernel(AdvectionRK4_3D) + pset.Kernel(zooplankton_drift)
pset.execute(kernels,   
             runtime=timedelta(hours=48),
#              endtime=datetime(2015,1,28,12,0,0),             
             dt=300,                       
             output_file=output_file)

output_file.close()
      