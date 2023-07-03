import fd3d
import numpy as np
import matplotlib.pyplot as plt

z50 = fd3d.run_sim()


half_range = np.arange(0, 1, 0.15)
contour_params = {
    'vmin': -1,
    'vmax': 1,
    'levels': np.concatenate((-half_range[1:][::-1], half_range))
}

for frame in [50, 100, 200, 300, 400, 500, 600, 700]:
    print(f'drawing frame {frame}')
    fig = plt.figure(figsize = (5,5))

    plt.xlabel('y')
    plt.ylabel('x')
    fig.subplots_adjust(
        bottom = 0.1, top = 0.9, 
        left = 0.1, right = 0.9
    )

    plt.contour(
        z50[frame // 2 - 1,:,:], **contour_params,
        colors = 'black', linewidths = 0.5, linestyles = 'solid'
    )
    plt.contourf(
        z50[frame // 2 - 1,:,:], **contour_params
    )

    # plt.colorbar(location = 'right', orientation = 'vertical', cax = ax)

    plt.vlines([40, 70], 0, 99, color = 'black')

    plt.savefig(f'img/snapshots/fd3d_z50_snap_{frame}')

    
