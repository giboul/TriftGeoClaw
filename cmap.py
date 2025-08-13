import matplotlib as mpl

water_color = (0.372, 0.588, 0.6)

def set_transparent_cmaps(eps):

    mpl.colormaps.register(mpl.colors.LinearSegmentedColormap("Reds_water", {
        'red': ((0.0, water_color[0], water_color[0]),
                (0.5, 1.0, 1.0,),
                (1.0, 1.0, 1.0,),),
        'green': ((0.0, water_color[1], water_color[1]),
                  (0.5, 1.0, 1.0,),
                  (1.0, 0.0, 0.0,),),
        'blue': ((0.0, water_color[2], water_color[2]),
                 (0.5, 1.0, 1.0,),
                 (1.0, 0.0, 0.0,),),
    }))


    mpl.colormaps.register(mpl.colors.LinearSegmentedColormap("RdBu_water", {
        'red': ((0.0, 1.0,1.0,),
                (0.2, 1.0, 1.0),
                (0.5, water_color[0], water_color[0]),
                (1.0, 0.0,0.0,),),
        'green': ((0.0, 0.0,0.0,), 
                (0.2, 0.0, 0.0),
                  (0.5, water_color[1], water_color[1]),
                  (1.0, 0.0,0.0,),),
        'blue': ((0.0, 0.0,0.0,),
                (0.2, 0.0, 0.0),
                 (0.5, water_color[2], water_color[2]),
                 (1.0, 1.0,1.0,),),
    }))

    mpl.colormaps.register(mpl.colors.LinearSegmentedColormap("RdBu_tc", {
        'red': ((0.0, 1.0, 1.0),
                (0.5, 1.0, 1.0),
                (1.0, 0.0, 0.0),),
        'green': ((0.0, 0.0, 0.0), 
                  (0.5, 1.0, 1.0),
                  (1.0, 0.0, 0.0),),
        'blue': ((0.0, 0.0, 0.0),
                 (0.5, 1.0, 1.0),
                 (1.0, 1.0, 1.0),),
        'alpha': ((0.0, 1.0, 1.0),
                  (0.5-eps, 1.0, 1.0),
                  (0.5, 0.0, 0.0),
                  (0.5+eps, 1.0, 1.0),
                  (1.0, 1.0, 1.0),),
    }))

    mpl.colormaps.register(mpl.colors.LinearSegmentedColormap("RdBu_tr", {
        'red': ((0.0, 1.0, 1.0),
                (0.5, 1.0, 1.0),
                (1.0, 0.0, 0.0),),
        'green': ((0.0, 0.0, 0.0), 
                  (0.5, 1.0, 1.0),
                  (1.0, 0.0, 0.0),),
        'blue': ((0.0, 0.0, 0.0),
                 (0.5, 1.0, 1.0),
                 (1.0, 1.0, 1.0),),
        'alpha': ((0.0, 0.0, 0.0),
                  (eps, 1.0, 1.0),
                  (0.5, 0.0, 0.0),
                  (1.0, 1.0, 1.0),),
    }))

    mpl.colormaps.register(mpl.colors.LinearSegmentedColormap("Blues_t", {
        'red': ((0.0, 1.0, 1.0),
                (1.0, 0.0, 0.0),),
        'green': ((0.0, 1.0, 1.0),
                  (1.0, 0.0, 0.0),),
        'blue': ((0.0, 1.0, 1.0),
                 (1.0, 1.0, 1.0),),
        'alpha': ((0.0, 0.0, 0.0),
                  (eps, 1.0, 1.0),
                  (1.0, 1.0, 1.0),),
    }))
