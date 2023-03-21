from numpy import *
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Cursor, LassoSelector, RadioButtons
import time


arr = array(range(6)) / 6.0
angle = 2 * pi * arr

cos_x = cos(angle)
sin_y = sin(angle)


class Ring:
    def __init__(self):
        self.x = cos_x
        self.y = sin_y
        self.z = zeros(6)
        self.phi = 180
        self.theta = 90
        self.q = 0.4
        self.center_mass = zeros(3)
        self.double_bond = 0
        self.bond_length = sqrt((self.x[0] - self.x[1]) ** 2 + (self.y[0] - self.y[1]) ** 2)
        self.oxygen = 2
        self.methil = 1
        self.double_methil = 5
        self.oh = 3
        self.type = "AST"
        self.plots = None

    def unmap(self):
        q2 = self.q * sin(pi * self.theta / 180)
        q3 = self.q * cos(pi * self.theta / 180)
        if self.type == "LUT":
            self.z = (1 / sqrt(3.0)) * q2 * cos(pi * self.phi / 180 + 2 * angle) + (1 / sqrt(6.0)) * q3 * (-1) ** array(range(6))
        elif self.type == "BCT":
            self.z = roll(flip(
                (1 / sqrt(3.0)) * q2 * cos(pi * self.phi / 180 + 2 * angle) + (1 / sqrt(6.0)) * q3 * (-1) ** array(
                    range(6)), 0), 1)
        else:
            self.z = flip((1 / sqrt(3.0)) * q2 * cos(pi * self.phi / 180 + 2 * angle) + (1 / sqrt(6.0)) * q3 * (-1) ** array(range(6)), 0)
        self.z = self.z - self.z[0]
        self.x = (0.5 + 0.5 * sqrt(1 - self.z ** 2)) * cos(angle)
        self.y = (0.5 + 0.5 * sqrt(1 - self.z ** 2)) * sin(angle)
        self.center_mass = array([
            self.x.sum() / 6,
            self.y.sum() / 6,
            self.z.sum() / 6,
        ])


ring = Ring()
ring.unmap()


def sinusoidal_projection(x, y):
    x = ((x - 180) * sin(pi * y / 180) + 180)
    return x, y


def update_figure():
    ring.unmap()
    [bonds, polyene_chain, polyene_db, methil, double_methil, double_bond, oxygen, oh, line] = ring.plots
    # cm.set_xdata([ring.center_mass[0]])
    # cm.set_ydata([ring.center_mass[1]])
    # cm.set_3d_properties([ring.center_mass[2]])
    # for at in range(6):
    #     atoms[at].set_xdata(ring.x)
    #     atoms[at].set_ydata(ring.y)
    #     atoms[at].set_3d_properties(ring.z)
    for bd in range(6):
        if bd < 5:
            bonds[bd].set_xdata([ring.x[bd], ring.x[bd + 1]])
            bonds[bd].set_ydata([ring.y[bd], ring.y[bd + 1]])
            bonds[bd].set_3d_properties([ring.z[bd], ring.z[bd + 1]])
        else:
            bonds[bd].set_xdata([ring.x[0], ring.x[bd]])
            bonds[bd].set_ydata([ring.y[0], ring.y[bd]])
            bonds[bd].set_3d_properties([ring.z[0], ring.z[bd]])

    double_bond.set_xdata([ring.x[ring.double_bond] * 0.8, ring.x[ring.double_bond + 1] * 0.8])
    double_bond.set_ydata([ring.y[ring.double_bond] * 0.8, ring.y[ring.double_bond + 1] * 0.8])
    double_bond.set_3d_properties([ring.z[ring.double_bond] * 0.8 - 0.04 * ring.bond_length, ring.z[ring.double_bond + 1] * 0.8 - 0.04 * ring.bond_length])

    if ring.type == "AST" or ring.type == "CAN":
        oxygen[0].set_xdata([ring.x[ring.oxygen] - 0.02 * ring.bond_length, ring.x[ring.oxygen] * 1.4 - 0.02 * ring.bond_length])
        oxygen[0].set_ydata([ring.y[ring.oxygen], ring.y[ring.oxygen] * 1.4])
        oxygen[0].set_3d_properties([ring.z[ring.oxygen] - 0.02 * ring.bond_length, ring.z[ring.oxygen] - 0.02 * ring.bond_length])

        oxygen[1].set_xdata([ring.x[ring.oxygen], ring.x[ring.oxygen] * 1.4])
        oxygen[1].set_ydata([ring.y[ring.oxygen] + 0.02 * ring.bond_length, ring.y[ring.oxygen] * 1.4 + 0.02 * ring.bond_length])
        oxygen[1].set_3d_properties([ring.z[ring.oxygen] + 0.02 * ring.bond_length, ring.z[ring.oxygen] + 0.02 * ring.bond_length])

        oxygen[2].set_xdata([ring.x[ring.oxygen] * 1.55])
        oxygen[2].set_ydata([ring.y[ring.oxygen] * 1.55])
        oxygen[2].set_3d_properties([ring.z[ring.oxygen]])

    if ring.type != "LUT":
        polyene_chain.set_xdata([ring.x[0], ring.x[0] + ring.bond_length, ring.x[0] + ring.bond_length + cos(ring.bond_length), ring.x[0] + ring.bond_length * 2 + cos(ring.bond_length)])
        polyene_chain.set_ydata([ring.y[0], 0, sin(ring.bond_length), ring.y[1]])
        polyene_chain.set_3d_properties([ring.z[0], ring.z[0], ring.z[0], ring.z[0]])
    else:
        polyene_chain.set_xdata(
            [ring.x[0], ring.x[0] + ring.bond_length, ring.x[0] + ring.bond_length + cos(ring.bond_length),
             ring.x[0] + ring.bond_length * 2 + cos(ring.bond_length)])
        polyene_chain.set_ydata([ring.y[0], 0, sin(ring.bond_length), ring.y[1]])
        polyene_chain.set_3d_properties([ring.z[0], ring.z[0] - ring.bond_length * 1, ring.z[0] - ring.bond_length * 2, ring.z[0] - ring.bond_length * 3])

    if ring.type != "LUT":
            polyene_db.set_xdata([ring.x[0] + ring.bond_length + ring.bond_length * 0.1 + sin(ring.bond_length * 0.13), ring.x[0] + ring.bond_length + cos(ring.bond_length) + ring.bond_length * 0.1]),
            polyene_db.set_ydata([0 - ring.bond_length * 0.1 + sin(ring.bond_length * 0.2), sin(ring.bond_length) - ring.bond_length * 0.1])
            polyene_db.set_3d_properties([ring.z[0], ring.z[0]]
    )
    else:
        polyene_db.set_xdata([ring.x[0] + ring.bond_length + ring.bond_length * 0.1 + sin(ring.bond_length * 0.13),
                              ring.x[0] + ring.bond_length + cos(ring.bond_length) + ring.bond_length * 0.1]),
        polyene_db.set_ydata(
            [0 - ring.bond_length * 0.1 + sin(ring.bond_length * 0.2), sin(ring.bond_length) - ring.bond_length * 0.1])
        polyene_db.set_3d_properties([ring.z[0] - ring.bond_length * 1.25, ring.z[0] - ring.bond_length * 2]
    )

    methil.set_xdata([ring.x[ring.methil], ring.x[ring.methil] * 1.6])
    methil.set_ydata([ring.y[ring.methil], ring.y[ring.methil] * 1.6])
    methil.set_3d_properties([ring.z[ring.methil], ring.z[ring.methil]])

    double_methil[0].set_xdata([ring.x[ring.double_methil], ring.x[ring.double_methil] * 1.4])
    double_methil[0].set_ydata([ring.y[ring.double_methil], ring.y[ring.double_methil] * 1.4])
    double_methil[0].set_3d_properties([ring.z[ring.double_methil], ring.z[ring.double_methil] + ring.bond_length * 0.4])

    double_methil[1].set_xdata([ring.x[ring.double_methil], ring.x[ring.double_methil] * 1.4])
    double_methil[1].set_ydata([ring.y[ring.double_methil], ring.y[ring.double_methil] * 1.4])
    double_methil[1].set_3d_properties([ring.z[ring.double_methil], ring.z[ring.double_methil] - ring.bond_length * 0.4])

    if ring.type == "AST" or ring.type == "LUT" or ring.type == "ZEA":
        oh[0].set_xdata([ring.x[ring.oh], ring.x[ring.oh] * 1.25])
        oh[0].set_ydata([ring.y[ring.oh], ring.y[ring.oh] * 1.25])
        oh[0].set_3d_properties([ring.z[ring.oh], ring.z[ring.oh] + ring.bond_length * 0.4 * (1 - 2 * (ring.type != "LUT"))])

        oh[1].set_xdata([ring.x[ring.oh] * 1.45])
        oh[1].set_ydata([ring.y[ring.oh] * 1.45])
        oh[1].set_3d_properties([ring.z[ring.oh] + ring.bond_length * 0.4 * (1 - 2 * (ring.type != "LUT"))])

    if ring.type == "AST" or ring.type == "CAN":
        line.set_xdata([ring.x[3], ring.x[5]])
        line.set_ydata([ring.y[3], ring.y[5]]),
        line.set_3d_properties([ring.z[3], ring.z[5]])
    elif ring.type == "LUT":
        line.set_xdata([ring.x[0], ring.x[3]])
        line.set_ydata([ring.y[0], ring.y[3]]),
        line.set_3d_properties([ring.z[0], ring.z[3]])
    else:
        line.set_xdata([ring.x[2], ring.x[5]])
        line.set_ydata([ring.y[2], ring.y[5]]),
        line.set_3d_properties([ring.z[2], ring.z[5]])

    img = plt.imread(f"{ring.type}.jpg")
    bg.set_data(img)

    pucker_map.set_xdata([sinusoidal_projection(ring.phi, ring.theta)[0]])
    pucker_map.set_ydata([sinusoidal_projection(ring.phi, ring.theta)[1]])
    fig.canvas.draw_idle()
    fig.canvas.flush_events()


# Define figure composition
fig = plt.figure(figsize=[10.67, 6])
# Define axes for 3d and map
ax_3d = fig.add_axes([0, 0.25, 0.5, 0.7], projection='3d')
ax_map = fig.add_axes([0.5, 0.25, 0.4, 0.7])
# Define axes for sliders
ax_phi = fig.add_axes([0.25, 0.15, 0.6, 0.03])
phi_slider = Slider(ax_phi, 'Phi', 0, 360, valinit=ring.phi)
ax_theta = fig.add_axes([0.25, 0.1, 0.6, 0.03])
theta_slider = Slider(ax_theta, 'Theta', 0, 180, valinit=ring.theta)
ax_q = fig.add_axes([0.25, 0.05, 0.6, 0.03])
q_slider = Slider(ax_q, 'Q', 0, 1, valinit=ring.q)

img = plt.imread(f"{ring.type}.jpg")
bg = ax_map.imshow(img, extent=[0, 360, 0, 180], aspect='auto')

rax = fig.add_axes([0.05, 0.05, 0.1, 0.13])
radio = RadioButtons(rax, ('AST', 'BCT', 'CAN', 'LUT', 'ZEA'))

# Create a map
[pucker_map] = ax_map.plot(
    *sinusoidal_projection(ring.phi, ring.theta), marker='o', color='red')

# Create 3d molecule
# [cm] = ax_3d.plot(*ring.center_mass, color='red', marker='o')
# atoms = []
# for atom in range(6):
#     atoms.append(
#         *ax_3d.plot(
#             [ring.x[atom]],
#             [ring.y[atom]],
#             [ring.z[atom]], color='black', marker=f'${atom}$')
#     )
def replot():
    bonds = []
    for bond in range(6):
        if bond < 5:
            bonds.append(*ax_3d.plot(
                [ring.x[bond], ring.x[bond + 1]],
                [ring.y[bond], ring.y[bond + 1]],
                [ring.z[bond], ring.z[bond + 1]], color='black', marker='', linewidth=2))
        else:
            bonds.append(*ax_3d.plot(
                [ring.x[0], ring.x[bond]],
                [ring.y[0], ring.y[bond]],
                [ring.z[0], ring.z[bond]], color='black', marker='', linewidth=2))

    if ring.type != "LUT":
        [polyene_chain] = ax_3d.plot(
            [ring.x[0], ring.x[0] + ring.bond_length, ring.x[0] + ring.bond_length + cos(ring.bond_length), ring.x[0] + ring.bond_length * 2 + cos(ring.bond_length)],
            [ring.y[0], 0, sin(ring.bond_length), ring.y[1]],
            [ring.z[0], ring.z[0], ring.z[0], ring.z[0]], color='black', marker='', linewidth=2
        )
    else:
        [polyene_chain] = ax_3d.plot(
            [ring.x[0], ring.x[0] + ring.bond_length, ring.x[0] + ring.bond_length + cos(ring.bond_length),
             ring.x[0] + ring.bond_length * 2 + cos(ring.bond_length)],
            [ring.y[0], 0, sin(ring.bond_length), ring.y[1]],
            [ring.z[0], ring.z[0] - ring.bond_length * 1, ring.z[0] - ring.bond_length * 2, ring.z[0] - ring.bond_length * 3], color='black', marker='', linewidth=2
        )
    if ring.type != "LUT":
        [polyene_db] = ax_3d.plot(
            [ring.x[0] + ring.bond_length + ring.bond_length * 0.1 + sin(ring.bond_length * 0.13), ring.x[0] + ring.bond_length + cos(ring.bond_length) + ring.bond_length * 0.1],
            [0 - ring.bond_length * 0.1 + sin(ring.bond_length * 0.2), sin(ring.bond_length) - ring.bond_length * 0.1],
            [ring.z[0], ring.z[0]], color='black', linewidth=2
    )
    else:
        [polyene_db] = ax_3d.plot(
            [ring.x[0] + ring.bond_length + ring.bond_length * 0.1 + sin(ring.bond_length * 0.13),
             ring.x[0] + ring.bond_length + cos(ring.bond_length) + ring.bond_length * 0.1],
            [0 - ring.bond_length * 0.1 + sin(ring.bond_length * 0.2), sin(ring.bond_length) - ring.bond_length * 0.1],
            [ring.z[0] - ring.bond_length * 1, ring.z[0] - ring.bond_length * 2], color='black', linewidth=2
    )

    [methil] = ax_3d.plot(
        [ring.x[ring.methil], ring.x[ring.methil] * 1.6],
        [ring.y[ring.methil], ring.y[ring.methil] * 1.6],
        [ring.z[ring.methil], ring.z[ring.methil]], color='black', linewidth=2
    )

    double_methil = [
        *ax_3d.plot(
            [ring.x[ring.double_methil], ring.x[ring.double_methil] * 1.4],
            [ring.y[ring.double_methil], ring.y[ring.double_methil] * 1.4],
            [ring.z[ring.double_methil], ring.z[ring.double_methil] + ring.bond_length * 0.4], color='black', linewidth=2),
        *ax_3d.plot(
            [ring.x[ring.double_methil], ring.x[ring.double_methil] * 1.4],
            [ring.y[ring.double_methil], ring.y[ring.double_methil] * 1.4],
            [ring.z[ring.double_methil], ring.z[ring.double_methil] - ring.bond_length * 0.4], color='black', linewidth=2),
    ]

# def replot_variables():
    if ring.type != "LUT":
        [double_bond] = ax_3d.plot(
            [ring.x[0] * 0.8, ring.x[1] * 0.8],
            [ring.y[0] * 0.8, ring.y[1] * 0.8],
            [ring.z[0] * 0.8 - 0.04 * ring.bond_length, ring.z[1] * 0.8 - 0.04 * ring.bond_length], color='black', linewidth=2
        )
    else:
        [double_bond] = ax_3d.plot(
            [ring.x[1] * 0.8, ring.x[2] * 0.8],
            [ring.y[1] * 0.8, ring.y[2] * 0.8],
            [ring.z[1] * 0.8 - 0.04 * ring.bond_length, ring.z[2] * 0.8 - 0.04 * ring.bond_length], color='black', linewidth=2
        )

    if ring.type == "AST" or ring.type == "CAN":
        oxygen = [
            *ax_3d.plot(
                [ring.x[ring.oxygen] - 0.02 * ring.bond_length, ring.x[ring.oxygen] * 1.4 - 0.02 * ring.bond_length],
                [ring.y[ring.oxygen], ring.y[ring.oxygen] * 1.4],
                [ring.z[ring.oxygen] - 0.02 * ring.bond_length, ring.z[ring.oxygen] - 0.02 * ring.bond_length], color='red', linewidth=2
            ),
            *ax_3d.plot(
                [ring.x[ring.oxygen], ring.x[ring.oxygen] * 1.4],
                [ring.y[ring.oxygen] + 0.02 * ring.bond_length, ring.y[ring.oxygen] * 1.4 + 0.02 * ring.bond_length],
                [ring.z[ring.oxygen] + 0.02 * ring.bond_length, ring.z[ring.oxygen] + 0.02 * ring.bond_length], color='red', linewidth=2
            ),
            *ax_3d.plot(
                [ring.x[ring.oxygen] * 1.55],
                [ring.y[ring.oxygen] * 1.55],
                [ring.z[ring.oxygen]], color='red', marker='$O$', linewidth=2, markersize=10
            )
        ]
    else:
        oxygen = None

    if ring.type == "AST" or ring.type == "LUT" or ring.type == "ZEA":
        oh = [
            *ax_3d.plot(
                [ring.x[ring.oh], ring.x[ring.oh] * 1.25],
                [ring.y[ring.oh], ring.y[ring.oh] * 1.25],
                [ring.z[ring.oh], ring.z[ring.oh] + ring.bond_length * 0.4 * (1 - 2 * (ring.type != "LUT"))], color='blue', linewidth=2),
            *ax_3d.plot(
                [ring.x[ring.oh] * 1.45],
                [ring.y[ring.oh] * 1.45],
                [ring.z[ring.oh] + ring.bond_length * 0.4 * (1 - 2 * (ring.type != "LUT"))], marker='$HO$', color='blue', markersize=20)
        ]
    else:
        oh = None

    if ring.type == "AST" or ring.type == "CAN":
        [line] = ax_3d.plot(
            [ring.x[3], ring.x[5]],
            [ring.y[3], ring.y[5]],
            [ring.z[3], ring.z[5]], color='black', linewidth=2, linestyle='--'
        )
    elif ring.type == "LUT":
        [line] = ax_3d.plot(
            [ring.x[0], ring.x[3]],
            [ring.y[0], ring.y[3]],
            [ring.z[0], ring.z[3]], color='black', linewidth=2, linestyle='--'
        )
    else:
        [line] = ax_3d.plot(
            [ring.x[2], ring.x[5]],
            [ring.y[2], ring.y[5]],
            [ring.z[2], ring.z[5]], color='black', linewidth=2, linestyle='--'
        )
    ring.plots = [bonds, polyene_chain, polyene_db, methil, double_methil, double_bond, oxygen, oh, line]


replot()

def onclick(event):
    # counter = 0
    # old_q = ring.q
    for point in event:
        # counter += 1
        ring.phi = ((point[0] - 180) / sin(pi * point[1] / 180) + 180)
        if ring.phi > 360:
            ring.phi = 360
        elif ring.phi < 0:
            ring.phi = 0
        ring.theta = point[1]
        if ring.theta > 180:
            ring.theta = 180
        elif ring.theta < 0:
            ring.theta = 0
        # if counter < len(event) / 2:
        #     ring.q -= old_q / (len(event) / 2)
        # elif counter > len(event) / 2 and counter != len(event) / 2:
        #     ring.q = old_q - ((len(event) - counter) * (old_q / (len(event) / 2)))
        # else:
        #     ring.q = old_q
        time.sleep(0.06)
        update_figure()


lasso = LassoSelector(ax=ax_map, useblit=False, onselect=onclick)


# Define figure updating
def sliders_on_changed(val):
    ring.phi = phi_slider.val
    ring.theta = theta_slider.val
    ring.q = q_slider.val
    update_figure()


def cartype(label):
    ring.type = label
    if ring.type != "LUT":
        ring.double_bond = 0
    else:
        ring.double_bond = 1
    ax_3d.clear()
    # double_bond.remove()
    # if oh is not None:
    #     oh[0].remove()
    #     oh[1].remove()
    # if oxygen is not None:
    #     oxygen[0].remove()
    #     oxygen[1].remove()
    #     oxygen[2].remove()
    # line.remove()
    replot()
    ax_3d.set_xlim(-1, 1)
    ax_3d.set_ylim(-1, 1)
    ax_3d.set_zlim(-1, 1)
    ax_3d.set_axis_off()
    update_figure()
    fig.canvas.draw()


radio.on_clicked(cartype)

phi_slider.on_changed(sliders_on_changed)
theta_slider.on_changed(sliders_on_changed)
q_slider.on_changed(sliders_on_changed)

# Map settings
ax_map.set_axis_off()
# ax_map.tick_params(axis='both', which='both', length=0)
# ax_map.spines['top'].set_visible(False)
# ax_map.spines['right'].set_visible(False)
# ax_map.spines['left'].set_position('center')
# ax_map.spines['bottom'].set_position('center')
# ax_map.spines['left'].set_visible(False)
# ax_map.spines['bottom'].set_visible(False)
# ax_map.xaxis.set_tick_params(bottom='on', top='off')
# ax_map.yaxis.set_tick_params(left='on', right='off')
# ax_map.plot([180, 180], [-1, 181], '-k')
# ax_map.plot([0, 360], [90, 90], '-k')
# ax_map.set_xticks([30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330])
# ax_map.set_yticks([15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165])
# ax_map.set_xlim([-5, 365])
# ax_map.set_ylim([-5, 185])
# y = linspace(0, 180, 720)
# ax_map.plot((-180) * sin(pi * y / 180) + 180, (y - 0.6) * 181.2 / 180, '-k', lw=4)
# ax_map.plot((+180) * sin(pi * y / 180) + 180, (y - 0.6) * 181.2 / 180, '-k', lw=4)
# for phi_t in [30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]:
#     ax_map.plot((phi_t - 180) * sin(pi * y / 180) + 180, y, '--k', lw=0.5)
# for theta_t in [15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165]:
#     ax_map.plot(linspace(-180 * sin(pi * theta_t / 180) + 180, +180 * sin(pi * theta_t / 180) + 180, 10), theta_t * ones(10), '--k', lw=0.5)

# Molecule 3d settings
# ax_3d.set_proj_type('ortho')
ax_3d.set_xlabel('X')
ax_3d.set_ylabel('Y')
ax_3d.set_zlabel('Z')
ax_3d.set_xlim(-1, 1)
ax_3d.set_ylim(-1, 1)
ax_3d.set_zlim(-1, 1)
ax_3d.set_axis_off()

fig.show()
plt.show()
