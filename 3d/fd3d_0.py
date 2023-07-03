import fd3d
import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera

z50 = fd3d.run_sim()

fig, ax = plt.subplots(figsize = (5,5))
ax.set_xlabel('y')
ax.set_ylabel('x')
# title = fig.text(0.5, 0.92, "")
# fig.tight_layout()
fig.subplots_adjust(
    bottom = 0.1, top = 0.9, 
    left = 0.1, right = 0.9
)

half_range = np.arange(0, 1, 0.15)
contour_params = {
    'vmin': -1,
    'vmax': 1,
    'levels': np.concatenate((-half_range[1:][::-1], half_range))
}

camera = Camera(fig)
for frame in range(z50.shape[0]):
    if frame % 10 == 0:
        print(f'drawing frame {frame}')

    # plt.subplots_adjust(
    #     bottom = 0.1, top = 0.9, 
    #     left = 0.1, right = 0.8
    # )
    # title.set_text(f'E_z(z = 50), step = {(frame + 1) * 2}')
    # fig.text(0.5, 0.92, f'E_z(z = 50), step = {(frame + 1) * 2}')
    # ax.set_title(f'E_z(z = 50), step = {(frame + 1) * 2}')

    ax.contour(
        z50[frame,:,:], **contour_params,
        colors = 'black', linewidths = 0.5, linestyles = 'solid'
    )
    ax.contourf(
        z50[frame,:,:], **contour_params
    )

    # plt.colorbar(location = 'right', orientation = 'vertical', cax = ax)

    ax.vlines([40, 70], 0, 99, color = 'black')
    
    camera.snap()

camera.animate().save(
    'img/fd3d_z50_0.gif', fps = 15, 
    writer = 'ffmpeg'
)
