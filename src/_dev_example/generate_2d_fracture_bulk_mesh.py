# %%
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt

import meshio
import pygmsh 
import gmsh

# %%
def drop_dupplicate_nodes(mesh, tol=1e-8):
    
    # Build a mapping old_id -> new_id after merging
    points = mesh.points
    N = len(points)
    new_index = -np.ones(N, dtype=int)

    # We’ll use a dict with rounded coords as keys
    coord_map = {}
    unique_points = []
    for i, p in enumerate(points):
        key = tuple(np.round(p / tol).astype(int))
        if key in coord_map:
            new_index[i] = coord_map[key]
        else:
            new_index[i] = len(unique_points)
            coord_map[key] = new_index[i]
            unique_points.append(p)

    unique_points = np.array(unique_points)

    print(f"Nodes dupplicates merging : {points.shape[0]=} vs {unique_points.shape[0]=}")

    new_cells = []
    for block in mesh.cells:
        new_data = np.vectorize(lambda idx: new_index[idx])(block.data)
        new_cells.append((block.type, new_data))
    
    return meshio.Mesh(
        points = unique_points,
        cells = new_cells,
        cell_data = {"gmsh:physical" : [np.arange(data.shape[0]) for _,data in new_cells]}
    )

# %%
def mesh_box_w_fracture(box_height = 1.0, box_width = 1.0, fracture_length = 0.5, box_res = 0.1, fracture_res=0.01):
    
    with pygmsh.occ.Geometry() as geom:
        
        # x0 = [-box_size/2, -box_size/2, 0.0]
        # bulk = geom.add_rectangle(x0, box_size, box_size)

        bulk_corner_top_right    = geom.add_point([ box_width/2, box_height/2, 0.0],  box_res)
        bulk_corner_top_left     = geom.add_point([-box_width/2, box_height/2, 0.0],  box_res)
        bulk_corner_bottom_right = geom.add_point([ box_width/2, -box_height/2, 0.0], box_res)
        bulk_corner_bottom_left  = geom.add_point([-box_width/2, -box_height/2, 0.0], box_res)
        
        box_edges = [
            geom.add_line(bulk_corner_top_left, bulk_corner_top_right),
            geom.add_line(bulk_corner_top_right, bulk_corner_bottom_right),
            geom.add_line(bulk_corner_bottom_right, bulk_corner_bottom_left),
            geom.add_line(bulk_corner_bottom_left, bulk_corner_top_left)
        ]
        
        bulk_loop = geom.add_curve_loop(box_edges)
        bulk_surface = geom.add_plane_surface(bulk_loop)
        
        frac_edge_left = geom.add_point([-fracture_length/2, 0.0, 0.0], fracture_res)
        frac_edge_right = geom.add_point([fracture_length/2, 0.0, 0.0], fracture_res)
        fracture = geom.add_line(frac_edge_left, frac_edge_right)
        
        # surfaces union
        new_surfaces = geom.boolean_fragments(bulk_surface, fracture)
        
        # Mark physical groups
        geom.add_physical(bulk_surface, label="bulk")
        geom.add_physical(fracture, label="fracture")
        
        geom.synchronize()
        
        geom.synchronize()
        mesh = geom.generate_mesh(order=1,algorithm=2)
        
    # # --- Extract only the fracture line elements
    # fracture_lines = mesh.cells_dict["line"]
    # fracture_line_data = mesh.cell_data_dict["gmsh:physical"]["line"]

    # # find the integer tag of "fracture"
    # fracture_tag = [tag for tag, name in zip(mesh.field_data.values(), mesh.field_data.keys()) if name == "fracture"][0][0]

    # # filter lines belonging to fracture
    # fracture_only = fracture_lines[fracture_line_data == fracture_tag]

    # Only keep triangles and fracture lines
    mesh =  meshio.Mesh(
        points = mesh.points,
        cells = {
            "line" : mesh.cells_dict['line'][mesh.cell_sets_dict["fracture"]["line"]],
            "triangle" : mesh.cells_dict['triangle']
            },
    )
        
    return drop_dupplicate_nodes(mesh, tol=1e-6)

# if __name__ == "__main__":

# Coarse mesh :

# %%
L_fracture = 100
L_bulk = L_fracture*1.5
res_bulk_outer = L_bulk/5.0
res_fracture = L_fracture/20.0

mesh = mesh_box_w_fracture(L_bulk/2, L_bulk, L_fracture, res_bulk_outer, res_fracture)

print(f"Coarse mesh : {mesh}")

# fig, ax = plt.subplots()

# # Plot the triangles
# coor = mesh.points
# conn = mesh.cells_dict['triangle']

# triang = matplotlib.tri.Triangulation(coor[:,0], coor[:,1], triangles=conn, mask=None)
# ax.triplot(triang, 'b-', lw=1)

# # Plot the lines
# conn = mesh.cells_dict['line']
# for i0, i1 in conn:
#     x = [coor[i0, 0], coor[i1, 0]]
#     y = [coor[i0, 1], coor[i1, 1]]
#     plt.plot(x, y, 'r-')  # black line

# ax.set_aspect('equal')

# plt.show()

#%% 

# Write to npy arrays

coor = mesh.points[:,:2]
tri_conn = mesh.get_cells_type("triangle")
seg_conn = mesh.cells_dict['line']

# Swap tri conn columns to make it CCW
tri_conn[:, [1, 2]] = tri_conn[:, [2, 1]]

np.save("mesh_triseg_coarse_0_coor", coor.astype(np.float64))
np.save("mesh_triseg_coarse_1_conn_tri", tri_conn.astype(np.int32))
np.save("mesh_triseg_coarse_2_conn_seg", seg_conn.astype(np.int32))

# Fine mesh

L_fracture = 100
L_bulk = L_fracture*1.5
res_bulk_outer = L_bulk/10.0
# res_fracture = L_fracture/1000.0
res_fracture = L_fracture/200.0

mesh = mesh_box_w_fracture(L_bulk/2, L_bulk, L_fracture, res_bulk_outer, res_fracture)

print(f"Fine mesh : {mesh}")

coor = mesh.points[:,:2]
tri_conn = mesh.get_cells_type("triangle")
seg_conn = mesh.cells_dict['line']

# Swap tri conn columns to make it CCW
tri_conn[:, [1, 2]] = tri_conn[:, [2, 1]]

np.save("mesh_triseg_fine_0_coor", coor.astype(np.float64))
np.save("mesh_triseg_fine_1_conn_tri", tri_conn.astype(np.int32))
np.save("mesh_triseg_fine_2_conn_seg", seg_conn.astype(np.int32))




