import math
import numpy as np

class Radar_prep(object):

    def __init__(self):
        print('hi')

    def ordered_az(self, var):

        if len(self.sweep_number['data']) == 11:
            ds = np.zeros(shape=(360*9, 1840),dtype=float)
        else:
            raise TypeError
        
        field_dict = {}
        if var in ['dBZ', 'lat', 'lon', 'alt']:
            level_loop = [0,2,4,5,6,7,8,9,10]
        elif var == 'VEL':
            level_loop = [1,3,4,5,6,7,8,9,10]
            
        for i_loop, isweep in enumerate(level_loop):
            sweep_slice = self.get_slice(isweep)
            hashTable = np.zeros((360,1840), dtype=float)  # Create an array with 360 slots, all initialized to 0
            for i_ray in range(sweep_slice.start, sweep_slice.stop):
                i_az=self.azimuth['data'][i_ray]
                if i_az >= 359.5:
                    i_az = i_az - 360
                if var == 'lat':
                    hashTable[round(i_az),:] = self.gate_latitude['data'][i_ray]
                elif var == 'lon':
                    hashTable[round(i_az),:] = self.gate_longitude['data'][i_ray]
                elif var == 'alt':
                    ele = [np.sin(math.radians(number)) for number in self.elevation['data'][i_ray].reshape(-1, 1) ]
                    hashTable[round(i_az),:] =  ele * self.range['data'].reshape(-1, 1).T
                else:
                    hashTable[round(i_az),:]= self.fields[var]['data'][i_ray]
                
            ds[i_loop*360:(i_loop+1)*360,:] = hashTable[:,:]
        field_dict['data'] = ds
        self.fields[var+'_reorder'] = field_dict
        return self

def height_exist(x, val:list, eps:int=50):
    '''
    to determine if the height has already been filled with obs.
    '''
    for i in val:
        if abs(x - i) <= eps:
            return True
    return False