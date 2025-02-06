import numpy as np
from numpy._core.defchararray import center
from vispy import scene, geometry
from vispy.scene.visuals import Mesh, Line, Sphere
from vispy.color import Color
from vispy.scene import visuals
from vispy.visuals.transforms import STTransform

# Parameters
detector_length = 200  # cm
detector_thickness = 2.54  # cm
detector_inner_radii = np.array([45] * 5) + 5 * np.array(range(5))  # MUST BE SORTED
detector_volume_inner_rad = 45  # cm
detector_volume_outer_rad = 75  # cm


def read_file(filename):
    with open(filename, "r") as file:
        file_content = file.read()
    groups = file_content.strip().split("\n\n")

    # Initialize the result lists
    result = []

    for group in groups[:3]:
        # Split each group into lines and convert to tuples of floats
        tuples = [tuple(map(float, line.split())) for line in group.strip().split("\n")]
        result.append(tuples)
    for group in groups[-2:]:
        detected = [int(bit) for bit in group if bit in "01"]
        result.append(detected)

    return result[0][0], result[1], result[2], result[3], result[4]


def create_cylinder_mesh(radius, length, num_segments):
    vertices = []
    faces = []
    edges = []

    # Generate vertices for the cylinder
    for i in range(num_segments):
        theta = 2 * np.pi * i / num_segments
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)

        # Bottom vertex
        vertices.append([x, y, -length / 2])
        # Top vertex
        vertices.append([x, y, length / 2])

    # Generate faces for the bottom and top caps
    bottom_center = len(vertices)
    top_center = bottom_center + 1
    vertices.append([0, 0, -length / 2])  # Bottom center
    vertices.append([0, 0, length / 2])  # Top center
    # Generate faces for the sides
    for i in range(num_segments):
        i0 = 2 * i
        i1 = 2 * ((i + 1) % num_segments)
        faces += [
            [i0, i1, i1 + 1],
            [i0, i1 + 1, i0 + 1],
            [bottom_center, i0, i1],
            [top_center, i1 + 1, i0 + 1],
        ]
        edges += [[i0 + 1, i1 + 1], [i0, i1], [i0, i0 + 1]]

    return np.array(vertices), np.array(faces), np.array(edges)


def draw_cylinder(radius, length, num_segments):
    # Create the cylinder mesh
    vertices, faces, edges = create_cylinder_mesh(radius, length, num_segments)

    # Create the solid cylinder mesh visual
    cylinder_mesh = Mesh(vertices=vertices, faces=faces, color=Color("lightblue"))
    wireframe = visuals.Line(
        pos=vertices, connect=edges, color=Color("black"), method="gl"
    )
    view.add(wireframe)
    view.add(cylinder_mesh)


# Create tube mesh
def create_tube_mesh(inner_radius, outer_radius, length, num_segments):
    vertices = []
    faces = []
    edges = []

    # Generate vertices for inner and outer surfaces
    for i in range(num_segments):
        theta = 2 * np.pi * i / num_segments
        x_inner = inner_radius * np.cos(theta)
        y_inner = inner_radius * np.sin(theta)
        x_outer = outer_radius * np.cos(theta)
        y_outer = outer_radius * np.sin(theta)

        # Bottom vertices
        vertices.append([x_inner, y_inner, -length / 2])
        vertices.append([x_outer, y_outer, -length / 2])
        # Top vertices
        vertices.append([x_inner, y_inner, length / 2])
        vertices.append([x_outer, y_outer, length / 2])

    # Generate faces for the sides
    for i in range(num_segments):
        i0 = 4 * i
        i1 = 4 * ((i + 1) % num_segments)
        # Faces
        faces += [
            [i0, i1, i1 + 2],
            [i0, i0 + 2, i1 + 2],
            [i0 + 1, i1 + 1, i1 + 3],
            [i0 + 1, i0 + 3, i1 + 3],
            [i0, i1, i0 + 1],
            [i1, i1 + 1, i0 + 1],
            [i0 + 2, i1 + 2, i0 + 3],
            [i1 + 2, i1 + 3, i0 + 3],
        ]
        # Edges
        edges += [
            [i0, i0 + 2],
            [i0 + 3, i0 + 1],
            [i0, i1],
            [i0 + 1, i1 + 1],
            [i0 + 2, i1 + 2],
            [i0 + 3, i1 + 3],
        ]

    return np.array(vertices), np.array(faces), np.array(edges)


# Create a canvas
# Create a canvas and view
canvas = scene.SceneCanvas(keys="interactive", bgcolor="white")
view = canvas.central_widget.add_view()
# Set up the camera
view.camera = scene.TurntableCamera(elevation=30, azimuth=30, distance=10000, fov=0.0)


def draw_tube(inner_radius, outer_radius, length, num_segments):
    # Create the tube mesh
    vertices, faces, edges = create_tube_mesh(
        inner_radius, outer_radius, length, num_segments
    )

    # Create the solid tube mesh visual
    tube_mesh = Mesh(vertices=vertices, faces=faces, color=Color("lightyellow"))
    wireframe = visuals.Line(
        pos=vertices, connect=edges, color=Color("black"), method="gl"
    )
    view.add(wireframe)
    view.add(tube_mesh)


def draw_sphere(center, color, radius=1):
    sphere = Sphere(
        radius,
        method="latitude",
        parent=view.scene,
        color=color,
    )
    sphere.transform = STTransform(translate=center)
    sphere.set_gl_state(
        "translucent",
        depth_test=False,
        blend=True,
        blend_func=("src_alpha", "one_minus_src_alpha"),
    )
    view.add(sphere)


def draw_path(points, detected):
    path = Line(pos=points, color="red", width=1)
    path.set_gl_state(
        "translucent",
        depth_test=False,
        blend=True,
        blend_func=("src_alpha", "one_minus_src_alpha"),
    )
    for i in range(1, len(points)):
        draw_sphere(points[i], Color("green" if detected[i] else "grey"))
    view.add(path)


def draw_annihilation(path1, detected1, path2, detected2, origin):
    draw_path([origin] + path1, [0] + detected1)
    draw_path([origin] + path2, [0] + detected2)
    draw_sphere(origin, Color("red"))


for inner_radius in detector_inner_radii:
    draw_tube(inner_radius, inner_radius + detector_thickness, detector_length, 30)
draw_cylinder(10.6, 4, 30)
origin, path1, path2, detected1, detected2 = read_file("visualization.data")
draw_annihilation(path1, detected1, path2, detected2, origin)
# Show the canvas
canvas.show()

# Run the application
if __name__ == "__main__":
    import sys

    if sys.flags.interactive != 1:
        canvas.app.run()
