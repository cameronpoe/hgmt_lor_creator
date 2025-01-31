import numpy as np
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


# Create cylinder mesh
def create_cylinder_mesh(inner_radius, outer_radius, height, num_segments):
    vertices = []
    faces = []

    # Generate vertices for inner and outer surfaces
    for i in range(num_segments):
        theta = 2 * np.pi * i / num_segments
        x_inner = inner_radius * np.cos(theta)
        y_inner = inner_radius * np.sin(theta)
        x_outer = outer_radius * np.cos(theta)
        y_outer = outer_radius * np.sin(theta)

        # Bottom vertices
        vertices.append([x_inner, y_inner, -height / 2])
        vertices.append([x_outer, y_outer, -height / 2])
        # Top vertices
        vertices.append([x_inner, y_inner, height / 2])
        vertices.append([x_outer, y_outer, height / 2])

    # Generate faces for the sides
    for i in range(num_segments):
        i0 = 4 * i
        i1 = 4 * ((i + 1) % num_segments)
        # Inner face
        faces.append([i0, i1, i0 + 2])
        faces.append([i1, i1 + 2, i0 + 2])
        # Outer face
        faces.append([i0 + 1, i1 + 1, i0 + 3])
        faces.append([i1 + 1, i1 + 3, i0 + 3])
        # Bottom face
        faces.append([i0, i1, i0 + 1])
        faces.append([i1, i1 + 1, i0 + 1])
        # Top face
        faces.append([i0 + 2, i1 + 2, i0 + 3])
        faces.append([i1 + 2, i1 + 3, i0 + 3])

    return np.array(vertices), np.array(faces)


# Create a canvas
# Create a canvas and view
canvas = scene.SceneCanvas(keys="interactive", bgcolor="white")
view = canvas.central_widget.add_view()
# Set up the camera
view.camera = scene.TurntableCamera(elevation=30, azimuth=30, distance=10000, fov=0.0)


def draw_cylinder(inner_radius, outer_radius, height, num_segments):
    # Create the cylinder mesh
    vertices, faces = create_cylinder_mesh(
        inner_radius, outer_radius, height, num_segments
    )

    # Create the solid cylinder mesh visual
    cylinder_mesh = Mesh(
        vertices=vertices, faces=faces, color=Color("lightblue", alpha=0.5)
    )
    # Make stuff see through
    view.add(cylinder_mesh)

    # Create wireframe lines along the length of the cylinder
    wireframe_vertices = []
    for i in range(num_segments):
        theta = 2 * np.pi * i / num_segments
        x_inner = inner_radius * np.cos(theta)
        y_inner = inner_radius * np.sin(theta)
        x_outer = outer_radius * np.cos(theta)
        y_outer = outer_radius * np.sin(theta)

        # Bottom to top lines for inner and outer surfaces
        wireframe_vertices.append([x_inner, y_inner, -height / 2])
        wireframe_vertices.append([x_inner, y_inner, height / 2])
        wireframe_vertices.append([x_outer, y_outer, -height / 2])
        wireframe_vertices.append([x_outer, y_outer, height / 2])

        # Horizontal lines at the bottom and top
        theta_next = 2 * np.pi * ((i + 1) % num_segments) / num_segments
        x_inner_next = inner_radius * np.cos(theta_next)
        y_inner_next = inner_radius * np.sin(theta_next)
        x_outer_next = outer_radius * np.cos(theta_next)
        y_outer_next = outer_radius * np.sin(theta_next)

        # Bottom inner and outer circles
        wireframe_vertices.append([x_inner, y_inner, -height / 2])
        wireframe_vertices.append([x_inner_next, y_inner_next, -height / 2])
        wireframe_vertices.append([x_outer, y_outer, -height / 2])
        wireframe_vertices.append([x_outer_next, y_outer_next, -height / 2])

        # Top inner and outer circles
        wireframe_vertices.append([x_inner, y_inner, height / 2])
        wireframe_vertices.append([x_inner_next, y_inner_next, height / 2])
        wireframe_vertices.append([x_outer, y_outer, height / 2])
        wireframe_vertices.append([x_outer_next, y_outer_next, height / 2])

    wireframe_vertices = np.array(wireframe_vertices)

    # Create the wireframe visual with thicker lines
    wireframe = Line(pos=wireframe_vertices, connect="segments", color="black", width=3)
    view.add(wireframe)


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
    draw_cylinder(inner_radius, inner_radius + detector_thickness, detector_length, 30)
origin = (-15.505946, 33.101066, 14.267038)
path1 = [(24.692972, 45.146507, 19.498425), (21.201054, 56.841789, 30.257233)]
path2 = [
    (-39.940784, 25.779266, 11.087149),
    (-39.570919, 22.475040, 9.805506),
    (52.030170, -6.920190, -28.633856),
    (51.959389, -6.434876, -28.646063),
    (55.376598, -5.862780, -32.032188),
    (60.668129, -0.803248, -37.574165),
]
detected1 = [1, 1]
detected2 = [1, 1, 1, 1, 1, 1]
draw_annihilation(path1, detected1, path2, detected2, origin)
# Show the canvas
canvas.show()

# Run the application
if __name__ == "__main__":
    import sys

    if sys.flags.interactive != 1:
        canvas.app.run()
