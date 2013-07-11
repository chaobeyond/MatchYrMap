##"""Library module for performing tasks useful for route modeling.
##
##To use this module, add the following lines to your code:
##import sys
##sys.path.append(r'T:\Models\BikeModels\LCOG\Scripts')
##from RouteModelingLib import *
##
##Note: Functions labelled with '_arcgis'require the local computer to have the
##library arcpy installed. This library is installed as part of Esri ArcGIS 10.x
##or greater.
##
##The best place to get the MySQLdb installer is at
##<codegood.com/archives/4>.
##
##Note: Functions labelled with '_mysql'require the local computer to have the
##library MySQLdb installed. The best place to get the MySQLdb installer is at
##<codegood.com/archives/4>.
##
##Note: Functions labelled with 'transform_' require the local computer to have
##the library pyproj installed. The best place to get the pyproj installer is at
##<code.google.com/p/pyproj>.
##
##Development code name "Breakin" (alphabetical).
##
##"""
##
### Loader info dictionaries.
##DEMO_LINK_NETWORK_FROM_ARCGIS = {
##    'source_path': 'feature_class_path_string',
##    'source_data_type': 'ArcGIS',  # Not yet used; just for info.
##    'source_id_field': 'id_field_string',
##    'source_from_node_field': 'from_node_field_name_string',
##    'source_to_node_field': 'to_node_field_name_string',
##    'source_from_level_field': 'from_level_field_name_string',
##    'source_to_level_field': 'to_level_field_name_string',
##    'shape_field': 'field_name_string',
##    'source_subselect_sql': 'sql_filter_string',
##    }
##
### Works (2012-07-06): Added link tree with walk-backs.
##class LinkNetwork():
##    """Stores information for a link network."""
##
##    def __init__(self, source_info={}):
##        self.source_path = source_info.get('source_path')
##        self.source_data_type = source_info.get('source_data_type')
##        self.source_id_field = source_info.get('source_id_field')
##        self.source_from_node_field = source_info.get('source_from_node_field')
##        self.source_to_node_field = source_info.get('source_to_node_field')
##        self.source_from_level_field = source_info.get('source_from_level_field')
##        self.source_to_level_field = source_info.get('source_to_level_field')
##        self.source_subselect_sql = source_info.get('source_subselect_sql')
##
##        self.link_info = {}
##        self.node_info = {}
##
##
##    # Links: Queries.
##
##    # Testing (2012-07-09): Deprecated max_node_distance.
##    def find_near_link_ids(self, link_id):
##        """Returns a list of all link IDs near a given link.
##
##        The given link_id is included in the result (that link is near itself).
##
##
##        """
##        nodes_to_cross = list(self.get_link_nodes(link_id))
##        near_link_ids = []
##
##        # Get the links at each node.
##        for node_id, node_level in nodes_to_cross:
##            node_links = self.get_link_ids_at_node(node_id, node_level)
##            # Add new link IDs to links_crossed.
##            for new_link_id in node_links:
##                if new_link_id not in near_link_ids:
##                    near_link_ids.append(new_link_id)
##
##        # Remove link ID from input.
##        try:
##            near_link_ids.remove(link_id)
##        except ValueError:
##            pass
##
##        return near_link_ids
##
##
##    # Works (2012-06-29): Deprecate.
##    def find_near_link_ids_old(self, link_id, max_node_distance=1):
##        """Returns a list of all link IDs near a given link.
##
##        The max_node_distance limits how many nodes can be crossed to reach a
##        near-link. If it is 1, only links that share a node with the given link
##        will be returned.
##
##        The given link_id is included in the result (that link is near itself).
##
##
##        """
##        crossed_nodes = 0
##        # Add values for the 0-link.
##        links_crossed = [link_id]
##        nodes_to_cross = list(self.get_link_nodes(link_id))
##        nodes_crossed = []
##
##        while crossed_nodes < max_node_distance:
##            # Get the links at each node.
##            for node_id, node_level in nodes_to_cross:
##                print('[find_near_link_ids] node = %i.' % node_id)
##                node_links = self.get_link_ids_at_node(node_id, node_level)
##                print('[find_near_link_ids] node_links = %s.' % str(node_links))
##
##                # Add new link IDs to links_crossed.
##                for new_link_id in node_links:
##                    print('[find_near_link_ids] new_link_id = %s.' % str(new_link_id))
##                    if new_link_id not in links_crossed:
##                        print 'not in links_crossed'
##                        links_crossed.append(new_link_id)
##                    for new_node in self.get_link_nodes(new_link_id):
##                        # Add node at other end of new link.
##                        if new_node not in nodes_to_cross:
##                            nodes_to_cross.append(new_node)
##                        else:
##                            nodes_to_cross.remove(new_node)
##
##            crossed_nodes += 1
##
##        # Remove link ID from input.
##        try:
##            links_crossed.remove(link_id)
##        except ValueError:
##            pass
##
##        return links_crossed
##
##
##    # Works (unknown date).
##    def get_link_count(self):
##        """Returns the number of links in the network."""
##
##        return len(link_info)
##
##
##    # Works (unknown date).
##    def get_link_ids(self):
##        """Returns a list of all link IDs in the network."""
##
##        return self.link_info.keys()
##
##
##    # Works (unknown date).
##    def get_link_info(self, link_id):
##        """Returns all link info for the given link ID."""
##
##        this_link_info = {}
##        for info in self.link_info[link_id]:
##            this_link_info[info] = self.link_info[link_id].get(info)
##
##        return this_link_info
##
##
##    # Works (2012-07-06): New code.
##    def get_near_link_tree(self, link_id, max_node_distance=1):
##        """Returns a tree-dictionary of all link IDs near a given link.
##
##        The tree-dictionary contains all the near link IDs as keys. The values
##        for these keys are a list with the walk-back path to the original link.
##        (1) The original link will be in the tree, but has an empty walk-back.
##        (2) A link adjacent to the original will have only the original link in
##        its walk-back.
##        (3) Further-out links will have walk-backs larger than one link. All
##        links in the walk-back will also be in the tree, and will always end
##        with the original link.
##        Example: {100: [], 101: [100], 102: [100], 201: [102, 100]}
##
##        The max_node_distance limits how many nodes can be crossed to reach a
##        near-link. If it is 1, only links that share a node with the given link
##        will be returned.
##
##        If max_node_distance = 0 (no nodes crossed), only the original link
##        will be included in the result.
##
##        """
##        ##print('[get_near_link_tree] Start.')
##        ##print('[get_near_link_tree] Link %i.' % link_id)
##
##        # Add value for the 0-link.
##        near_link_tree = {link_id: []}
##        ##print('[get_near_link_tree] Added 0-link %i.' % link_id)
##        # Set up node collectors.
##        crossed_count = 0
##        nodes_to_cross = list(self.get_link_nodes(link_id))
##        nodes_crossed = []
##        ##print('[get_near_link_tree] Nodes to cross %s.' % str(nodes_to_cross))
##
##        while crossed_count < max_node_distance:
##            # Evaluate links at each node.
##            for node_id, node_level in  sorted(nodes_to_cross):
##                link_ids = self.get_link_ids_at_node(node_id, node_level)
##                for link_id in link_ids:
##                    # Add new links if ID not a key in tree.
##                    if link_id not in near_link_tree:
##                        ##print('[get_near_link_tree] Link %i not in tree.' % link_id)
##                        links_near = self.find_near_link_ids(link_id)
##                        ##print link_id, 'links near', links_near
##                        for near_link_id in links_near:
##                            ##print 'near link', near_link_id
##                            # Append the ID & it's tree-chain to the new link.
##                            if near_link_id in near_link_tree:
##                                near_link_tree[link_id] = (
##                                    [near_link_id] + near_link_tree[near_link_id]
##                                    )
##                                ##print('[get_near_link_tree] Added %i branch %s.' % (link_id, str(near_link_tree[link_id])))
##                                break
##                    for new_node in self.get_link_nodes(link_id):
##                        # Add node at other end of new link.
##                        if new_node not in nodes_to_cross and nodes_crossed:
##                            nodes_to_cross.append(new_node)
##                            ##print('[get_near_link_tree] Added node %s to nodes_to_cross.' % str(new_node))
##                nodes_crossed.append((node_id, node_level))
##                ##print('[get_near_link_tree] Removed node %i at level %i.' % (node_id, node_level))
##            crossed_count += 1
##            ##print('[get_near_link_tree] Crossed count = %i.' % crossed_count)
##
##        return near_link_tree
##
##
##    # Links: Data management.
##
##    # Works (unknown date).
##    def add_link_to_network(self, link_id, from_node, to_node,
##                            from_level=None, to_level=None, xyz_list=None):
##        """Adds a link to the network with the provided information."""
##
##        # Add track ID as key if not there yet.
##        self.link_info[link_id] = {
##            'from_node': from_node, 'to_node': to_node,
##            'from_level': from_level, 'to_level': to_level,
##            'points': {}
##                }
##        position = 1
##        if xyz_list:
##            for xyz in xyz_list:
##                self.link_info[link_id]['points'][position] = {
##                    'x': xyz[0], 'y': xyz[1], 'z': xyz[2]
##                    }
##                position += 1
##        self.add_node_to_network(from_node, from_level, [link_id])
##        self.add_node_to_network(to_node, to_level, [link_id])
##
##        return
##
##
##    # Works (unknown date).
##    def add_vertex_to_link(self, link_id, xyz, position=None):
##        """Adds provided vertex info to a link in the network.
##
##        If position is not noted, will be added at the end.
##
##        """
##        if position is None:
##            position = len(self.link_info[link_id]['points']) + 2
##        self.link_info[link_id]['points'][position] = {
##            'x': xyz[0], 'y': xyz[1], 'z': xyz[2]
##            }
##
##        return
##
##
##    # Nodes: Queries.
##
##    # Works (2012-06-29): Added logic to handle incorrect ID or level.
##    def get_link_ids_at_node(self, node_id, node_level=None):
##        """Returns a list of link IDs for a given node & level.
##
##        If no level is provided, links IDs at all levels of that node will be
##        returned.
##
##        """
##        # At the least, return an empty list if node_id not there, or
##        # node_level not on node.
##        link_ids = []
##
##        if node_id in self.node_info:
##            if node_level:
##                link_ids = self.node_info[node_id].get(node_level)
##            else:
##                for level in self.node_info[node_id]:
##                    link_ids.extend(self.node_info[node_id].get(level))
##
##        return link_ids
##
##
##    # Works (2012-06-29): Removed list conversion for result (now tuple).
##    def get_link_nodes(self, link_id):
##        """Returns a pair of tuples containing node IDs & levels."""
##
##        if link_id is None:
##            result = ((None, None), (None, None))
##
##        else:
##            fnode = (self.link_info[link_id]['from_node'],
##                     self.link_info[link_id]['from_level'])
##            tnode = (self.link_info[link_id]['to_node'],
##                     self.link_info[link_id]['to_level'])
##            result = fnode, tnode
##
##        return result
##
##
##    # Works (2012-06-29).
##    def get_node_xyz(self, node_id):
##        """Returns a tuple of the XYZ location for a given node."""
##
##        if node_id is None:
##            x, y , z = None, None, None
##        else:
##
##            node_links = self.get_link_ids_at_node(node_id)
##
##            a_link = node_links[0]
##
##            if node_id == self.link_info[a_link]['from_node']:
##                link_position = 1
##            elif node_id == self.link_info[a_link]['to_node']:
##                link_position = sorted(self.link_info[a_link]['points'])[-1]
##            # Recursively impossible, but we'll include it.
##            else:
##                x, y, z = None, None, None
##
##            x = self.link_info[a_link]['points'][link_position]['x']
##            y = self.link_info[a_link]['points'][link_position]['y']
##            z = self.link_info[a_link]['points'][link_position]['z']
##
##        return x, y, z
##
##
##    # Nodes: Data management.
##
##    # Works (unknown date).
##    def add_node_to_network(self, node_id, node_level=None, link_ids=[]):
##        """Add node ID as key if not there yet. Add level key if not there yet.
##
##        """
##        if node_id not in self.node_info:
##            self.node_info[node_id] = {}
##        if node_level not in self.node_info[node_id]:
##            self.node_info[node_id][node_level] = []
##        self.node_info[node_id][node_level].extend(link_ids)
##
##        return
##
##
##    # Vertices: Queries.
##
##    # Works (2012-05-31): New function.
##    def get_link_vertex(self, link_id, vertex_position):
##        """Returns an xyz tuple for a vertex on a link."""
##
##        # Test for validity of link ID & vertex position.
##        if link_id is None:
##            link_vertex = (None, None, None)
##        else:
##            if self.link_info[link_id]['points'].get(vertex_position) is None:
##                link_vertex = (None, None, None)
##
##            # All valid, assign tuple data.
##            else:
##                link_vertex = (
##                    self.link_info[link_id]['points'][vertex_position].get('x'),
##                    self.link_info[link_id]['points'][vertex_position].get('y'),
##                    self.link_info[link_id]['points'][vertex_position].get('z')
##                    )
##
##        return link_vertex
##
##
##    # Works (2012-05-31): Changed to use get_link_vertex().
##    def get_link_vertex_tuples(self, link_ids):
##        """Returns all the given links' vertices as a list of tuples.
##
##        Format: [(link_id, vertex_position, x, y, z), ...]
##
##        """
##        vertex_tuples = []
##        for link_id in link_ids:
##            for vertex_position in sorted(self.link_info[link_id]['points']):
##                vertex_tuples.append(
##                    (link_id, vertex_position,
##                     self.get_link_vertex(link_id, vertex_position))
##                    )
##
##        return vertex_tuples
##
##
##    # Works (2012-05-31): Changed to use get_link_vertex().
##    def get_link_vertices(self, link_id):
##        """Returns a list of xyz tuples for all vertices on a link."""
##
##        link_vertices = []
##
##        for vertex_position in sorted(self.link_info[link_id]['points']):
##            link_vertices.append(
##                self.get_link_vertex(link_id, vertex_position)
##                )
##
##        return link_vertices
##
##
### Works (2012-05-07).
##def load_link_network_from_arcgis(data_dictionary):
##    """Loads the link network from an ArcGIS-native spatial data type.
##
##    Requires a conversion dictionary for the source data. Returns a dictionary
##    containing each individual link as a sub-dictionary.
##
##    """
##    try:
##        import arcpy
##    except ImportError:
##        raise Exception
##
##    network = LinkNetwork(data_dictionary)
##
##    try:
##        ncursor = arcpy.SearchCursor(network.source_path,
##                                     network.source_subselect_sql)
##    except arcpy.ExecuteError:
##        print 'Need to add a error prompt here.'
##
##    for n in ncursor:
##        if network.source_from_level_field:
##            fzlev = n.getValue(network.source_from_level_field)
##        if network.source_to_level_field:
##            tzlev = n.getValue(network.source_to_level_field)
##
##        network.add_link_to_network(
##            link_id = n.getValue(network.source_id_field),
##            from_node = n.getValue(network.source_from_node_field),
##            to_node = n.getValue(network.source_to_node_field),
##            from_level = fzlev,
##            to_level = tzlev,
##            xyz_list = None
##            )
##
##        # Get the points from the link's geometry (only set for single-part).
##        point_position = 1
##        for point in n.getValue(data_dictionary['shape_field']).getPart(0):
##            network.add_vertex_to_link(
##                link_id = n.getValue(network.source_id_field),
##                xyz = (point.X, point.Y, point.Z), position = point_position
##                )
##            point_position += 1
##
##    del ncursor
##
##    return network