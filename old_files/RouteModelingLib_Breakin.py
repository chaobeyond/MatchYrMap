"""Library module for performing tasks useful for route modeling.

To use this module, add the following lines to your code:
import sys
sys.path.append(r'T:\Models\BikeModels\LCOG\Scripts')
from RouteModelingLib import *

Note: Functions labelled with '_arcgis'require the local computer to have the
library arcpy installed. This library is installed as part of Esri ArcGIS 10.x
or greater.

The best place to get the MySQLdb installer is at
<codegood.com/archives/4>.

Note: Functions labelled with '_mysql'require the local computer to have the
library MySQLdb installed. The best place to get the MySQLdb installer is at
<codegood.com/archives/4>.

Note: Functions labelled with 'transform_' require the local computer to have
the library pyproj installed. The best place to get the pyproj installer is at
<code.google.com/p/pyproj>.

Development code name "Breakin" (alphabetical).

"""

import datetime
from exceptions import Exception

# Loader info dictionaries.
DEMO_GPS_TRIP_FROM_MYSQL = {
    'source_host': 'database.domain-name-string.org',
    'source_port': 3306,  # Port as integer.
    'source_user': 'username_string',
    'source_password': 'user_password_string',
    'source_db': 'database_name_string',
    'source_table': 'table_name_as_string',
    'source_data_type': 'MySQL',  # Not yet used; just for info.
    'source_id_field': 'id_field_string',
    'source_timestamp_field': 'timestamp_field_name_string',
    'source_easting_field': 'easting_field_name_string',
    'source_northing_field': 'northing_field_name_string',
    'source_elevation_field': 'elevation_field_name_string',
    'source_xy_accuracy_field': 'xy_accuracy_field_name_string',
    'source_z_accuracy_field': 'z_accuracy_field_name_string',
    }
DEMO_LINK_NETWORK_FROM_ARCGIS = {
    'source_path': 'feature_class_path_string',
    'source_data_type': 'ArcGIS',  # Not yet used; just for info.
    'source_id_field': 'id_field_string',
    'source_from_node_field': 'from_node_field_name_string',
    'source_to_node_field': 'to_node_field_name_string',
    'source_from_level_field': 'from_level_field_name_string',
    'source_to_level_field': 'to_level_field_name_string',
    'shape_field': 'field_name_string',
    'source_subselect_sql': 'sql_filter_string',
    }
# Output updater data dictionaries.
DEMO_OD_PAIR_UPDATE_TO_ARCGIS = {
    'output_path': 'output_path_string',
    'trip_id_field': ('field_name_string', 'field_type_sting'),
    'node_id_field': ('field_name_string', 'field_type_sting'),
    'node_type_field': ('field_name_string', 'field_type_sting'),
    'shape_field': 'field_name_string'
    }
DEMO_OD_PAIR_UPDATE_TO_CSV = {
    'output_path': 'output_path_string',
    'trip_id_field': ('field_name_string', 'field_type_sting'),
    'node_id_field': ('field_name_string', 'field_type_sting'),
    'node_type_field': ('field_name_string', 'field_type_sting'),
    }
DEMO_ROUTES_UPDATE_TO_ARCGIS = {
    'output_path': 'output_path_string',
    'route_id_field': ('field_name_string', 'field_type_sting'),
    'link_id_field': ('field_name_string', 'field_type_sting'),
    'link_position_field': ('field_name_string', 'field_type_sting'),
    'shape_field': 'field_name_string',
    }
DEMO_ROUTES_UPDATE_TO_CSV = {
    'output_path': 'output_path_string',
    'route_id_field': ('field_name_string', 'field_type_sting'),
    'link_id_field': ('field_name_string', 'field_type_sting'),
    'link_position_field': ('field_name_string', 'field_type_sting')
    }
DEMO_ROUTEPOINTS_UPDATE_TO_ARCGIS = {
    'output_path': 'output_path_string',
    'route_id_field': ('field_name_string', 'field_type_sting'),
    'routepoint_id_field': ('field_name_string', 'field_type_sting'),
    # Optional.
    'link_ids_field': ('field_name_string', 'field_type_sting'),
    # Optional.
    'link_vertex_field': ('field_name_string', 'field_type_sting'),
    'shape_field': 'field_name_string'
    }
DEMO_TRACKPOINTS_UPDATE_TO_ARCGIS = {
    'output_path': 'output_path_string',
    'track_id_field': ('field_name_string', 'field_type_sting'),
    'trackpoint_id_field': ('field_name_string', 'field_type_sting'),
    # Optional.
    'link_id_field': ('field_name_string', 'field_type_sting'),
    # Optional.
    'link_vertex_field': ('field_name_string', 'field_type_sting'),
    'shape_field': 'field_name_string'
    }
DEMO_OUTPUT_LAYERS_CREATE_TO_ARCGIS = {
    'folder_path': 'folder_path_string',
    'name_format': 'layer_<ID>_file_name_string',  # No .lyr; added in script.
    'template_path': 'layer_path_string',
    'sublayers_sql': {
        'layer_name_string': 'definition_query_sql_string',  # Use <ID> for sub.
        'layer_name_string': 'definition_query_sql_string',  # Use <ID> for sub.
        'layer_name_string': 'definition_query_sql_string',  # Use <ID> for sub.
        'layer_name_string': 'definition_query_sql_string',  # Use <ID> for sub.
        }
    }


# Works (2012-07-06): Added multi-link routepoints; added debug notes.
class GPSTrip():
    """Stores trackpoint & other information for a GPS-recorded trip."""


    # Works (2012-05-21).
    def __init__(self, source_info={}):

        self.trip_id = source_info.get('trip_id')

        self.source_path = source_info.get('source_path')
        self.source_data_type = source_info.get('source_data_type')
        self.source_id_field = source_info.get('source_id_field')
        self.source_timestamp_field = source_info.get('source_timestamp_field')
        self.source_elevation_field = source_info.get('source_elevation_field')
        self.source_xy_accuracy_field = source_info.get('source_xy_accuracy_field')
        self.source_z_accuracy_field = source_info.get('source_z_accuracy_field')

        self.statistics = {}

        self.trackpoints = {}
        self.track_speeds = {}
        self.routepoints = {}
        self.routelinks = {}

        self.all_notes = []  # Repository for debugging notes.
        self.echo = False  # Will print debugging notes to the console.
        self.note('[__init__] trip_id %i.' % self.trip_id)

    def note(self, note_to_add):
        """Adds a note to self.all_notes."""

        timestamp = str(datetime.datetime.now().replace(microsecond = 0))

        self.all_notes.append('%s: %s' % (timestamp, note_to_add))

        if self.echo:
            print note_to_add

        return


    # Tracks: Queries.

    # Works (2012-07-09).
    def find_absolute_trip_length(self):
        """Calculates the absolute length of the trip.

        Absolute distance is a pure calculation of distance traveled between
        the ordered trackpoints in the trip. This is *not* the distance
        along the route as assigned by the map-matching.

        Returns a tuple containing:
        (absolute_distance, absolute_time_seconds)

        """

        # Init travel tracker.
        distance_traveled = 0.0
        time_traveled_seconds = 0.0

        for this_id in self.trackpoints:

            if this_id != 1:
                that_id = this_id - 1

                # This trackpoint's geometry.
                this_x = self.trackpoints[this_id].get('x')
                this_y = self.trackpoints[this_id].get('y')
                this_time = self.trackpoints[this_id].get('timestamp')

                # Previous trackpoint's geometry.
                that_x = self.trackpoints[that_id].get('x')
                that_y = self.trackpoints[that_id].get('y')
                that_time = self.trackpoints[that_id].get('timestamp')

                # Find distance & time between points.
                this_distance = (
                    ((that_x - this_x) ** 2 + (that_y - this_y) ** 2) ** 0.5
                    )
                this_time_seconds = (this_time - that_time).seconds

                # Add distance & time to total.
                distance_traveled += this_distance
                time_traveled_seconds += this_time_seconds

        return distance_traveled, time_traveled_seconds


    # Works (2012-05-31).
    def find_largest_trackpoint_gap(self):
        """Determines the largest gap between adjacent trackpoints.

        Returns a tuple containing:
        (largest_gap_distance, largest_gap_time_seconds)

        """

        # Init gap tracker.
        largest_gap_distance = 0.0
        largest_gap_time_seconds = 0.0

        for this_id in self.trackpoints:

            if this_id != 1:
                that_id = this_id - 1

                # This trackpoint's geometry.
                this_x = self.trackpoints[this_id].get('x')
                this_y = self.trackpoints[this_id].get('y')
                this_time = self.trackpoints[this_id].get('timestamp')

                # Previous trackpoint's geometry.
                that_x = self.trackpoints[that_id].get('x')
                that_y = self.trackpoints[that_id].get('y')
                that_time = self.trackpoints[that_id].get('timestamp')

                this_gap_distance = (
                    ((that_x - this_x) ** 2 + (that_y - this_y) ** 2) ** 0.5
                    )

                # If gap_distance largest sor far, assign.
                if this_gap_distance > largest_gap_distance:
                    largest_gap_distance = this_gap_distance
                    largest_gap_time_seconds = (this_time - that_time).seconds

        return largest_gap_distance, largest_gap_time_seconds


    # Doesn't currently work (2012-07-09).
    def find_trackpoints_outside_area(self, easting_range, northing_range):
        """Compiles all trackpoints that fall outside the given area.

        The easting & northing ranges represent an envelope that defines the
        "area".

        Returns a list of trackpoint IDs that are outside the envelope.

        """

        # Init ID list.
        trackpoints_outside = []

        for trackpoint_id in self.trackpoints:

            # This trackpoint's geometry.
            x = self.trackpoints[trackpoint_id].get('x')
            y = self.trackpoints[trackpoint_id].get('y')

            x_is_outside = x < easting_range[0] or x > easting_range[1]
            y_is_outside = y < northing_range[0] or y > northing_range[1]

            if x_is_outside or y_is_outside:
                trackpoints_outside.append(trackpoint_id)

        return trackpoints_outside


    # Works (2012-05-21).
    def get_trackpoint_count(self):
        """Returns the number of trackpoint in the trip."""

        return len(self.trackpoints)


    # Works (2012-05-25).
    def get_trackpoint_ids(self):
        """Returns all trackpoint IDs."""

        trackpoint_id_list = []

        for position in self.trackpoints:
            trackpoint_id_list.append(position)

        return trackpoint_id_list


    # Works (2012-05-25).
    def get_trackpoint_at(self, position):
        """Returns a trackpoint's info."""

        trackpoint_info = {}

        for key in self.trackpoints.get(position):
            trackpoint_info[key] = self.trackpoints[position].get(key)

        return trackpoint_info


    # Works (2012-05-25).
    def get_trackpoints(self):
        """Returns all trackpoint info."""

        trackpoint_list = []

        for position in self.trackpoints:
            trackpoint_list.append(self.trackpoints.get(position))

        return trackpoint_list


    # Tracks: Data management.

    # Works (2012-05-21).
    def add_trackpoint(self, position, timestamp, xyz, accuracy):
        """Adds a single trackpoint to the trip's dictionary"""

        self.trackpoints[position] = {
            'trackpoint_id': position, 'timestamp': timestamp,
            'x': xyz[0], 'y': xyz[1], 'z': xyz[2],
            'xy_accuracy': accuracy[0], 'z_accuracy': accuracy[1]
            }

        return position, xyz


    # Not sure this works at all (2012-05-21).
    def identify_problem_trackpoints(self):
        """Identifies & marks suspicious trackpoints via specific QA methods."""

        speed_mean = self.statistics['speed_mean']
        speed_stdev = self.statistics['speed_standard_deviation']

        for position in sorted(self.trackpoints):
            # Initial position is considered gospel. If bad, need to ID before
            # running most of these class functions (and map-matching).
            if position > 1:
                ##lastposition = position - 1
                thispoint = self.trackpoints[position]
                ##lastpoint = self.trackpoints[lastposition]
                # Move another point back to until no previous problems, or
                # you've reached the initial point (which is gospel).
                ##while lastpoint.get('problems') and lastposition != 1:
                ##    lastposition = lastposition - 1
                ##    lastpoint = self.trackpoints[lastposition]

                # Problem: Unlikely speed.
                # If speed is more than two standard deviations from the mean,
                # mark.
                deviation = abs(speed_mean - self.track_speeds['speed'])
                if deviation >= speed_stdev * 2:
                    if thispoint['problems']:
                        thispoint['problems'].append('speed')
                    else:
                        thispoint['problems'] = ['speed']

        return


    # Not yet written (2012-05-21).
    def reset_trackpoint_positions(self):
        """Resets trackpoint position values based on timestamp values.

        Still needs to be written.

        """

        return


    # Not working, unused (2012-005-21).
    def set_trackpoint_speed(self):
        """Assigns a speed to a trackpoint based on the previous one."""

        for position in self.trackpoints:
         # Initial position gets a speed of 0 (no previous point, plus logic
         # says you aren't moving until you've begun).
         if position == 1:
            self.track_speeds[position] = 0.0
         else:
            this_time = self.trackpoints[position]['timestamp']
            that_time = self.trackpoints[position - 1]['timestamp']
            dt = this_time - that_time
            time_elapsed = dt.seconds + dt.microseconds * 10 ** -6
            x2 = self.trackpoints[position]['x']
            y2 = self.trackpoints[position]['y']
            x1 = self.trackpoints[position - 1]['x']
            y1 = self.trackpoints[position - 1]['y']
            distance_elapsed = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5

            self.track_speeds[position] = (distance_elapsed / time_elapsed)

        # Update speed statistics.

        # Mean speed.
        speed_sum = 0
        speed_count = 0
        for position in self.trackpoints:
            # Skip the initial point for speed stats; logic precludes this.
            if position != 1:
                speed_sum += self.track_speeds[position]
                speed_count += 1
        speed_mean = speed_sum / speed_count

        # Standard deviation.
        for position in self.trackpoints:
            # Skip the initial point for speed stats; logic precludes this.
            variance_numer = 0
            variance_denom = 0
            if position != 1:
                variance_numer += (
                    self.track_speeds[position] - speed_mean) ** 2
                variance_denom += 1
        speed_variance = variance_numer / variance_denom
        speed_stdev = variance ** 0.5

        self.statistics['speed_mean'] = speed_mean
        self.statistics['speed_standard_deviation'] = speed_stdev

        return speed_mean


    # Routes: Queries.

    # Works (2012-07-06): Added multi-link routepoints.
    def find_routepoint(self, trackpoint_id, link_network, cursor_size=1):
        """Finds the nearest vertex in the link network.

        Previous routepoint must have been already solved.

        """
        self.note('[find_routepoint] Trackpoint %i.' % trackpoint_id)
        self.note('[find_routepoint] Info %s.' % str(self.get_trackpoint_at(trackpoint_id)))

        # Update potential link IDs.

        # If initial trackpoint, all links are fair game. Same goes if
        # something else weird happens (in that case, link_ids will be [None]).
        if trackpoint_id == 1:
            self.note('[find_routepoint] Is first trackpoint.')
            candidate_links = link_network.get_link_ids()
            candidate_links_tree = {}
            for candidate_link_id in candidate_links:
                candidate_links_tree[candidate_link_id] = []
        # Otherwise, get only links sharing a node with the cursor link.
        else:
            last_routepoint_id = trackpoint_id - 1
            last_link_ids = self.routepoints[last_routepoint_id].get('link_ids', [])
            while len(last_link_ids) == 0:
                # You're back at the startpoint and have no link_id. All links.
                if last_routepoint_id == 0:
                    self.note('[find_routepoint] Back to first routepoint.')
                    candidate_links = link_network.get_link_ids()
                    candidate_links_tree = {}
                    for candidate_link_id in candidate_links:
                        candidate_links_tree[candidate_link_id] = []
                    break
                # Move back another spot.
                else:
                    self.note('[find_routepoint] Last routepoint has no link ID: backing up.')
                    last_routepoint_id -= 1
                    last_link_ids = self.routepoints[last_routepoint_id].get('link_ids', [])[:1]
            # Contra argument would be if you've got no previous routepoint. (?)
            if len(last_link_ids) > 0:
                self.note('[find_routepoint] Last routepoint has link %i.' % last_link_ids[0])
                candidate_links_tree = link_network.get_near_link_tree(
                    link_id = last_link_ids[0],
                    max_node_distance = cursor_size
                    )
                candidate_links = candidate_links_tree.keys()

        self.note('[find_routepoint] Candidate links %s.' % str(candidate_links))
        self.note('[find_routepoint] Candidate links tree %s.' % str(candidate_links_tree))

        # Update the potential link vertices.
        candidate_vertices = link_network.get_link_vertex_tuples(candidate_links)

        # Find nearest candidate vertex.

        # Trackpoint geometry.
        xtp = self.trackpoints[trackpoint_id].get('x')
        ytp = self.trackpoints[trackpoint_id].get('y')
        ztp = self.trackpoints[trackpoint_id].get('z')

        # Init near evaluators.
        near_link_id = None
        near_vertex_position = None
        near_xyz = None
        near_distance = None

        for vertex in candidate_vertices:
            xcv, ycv, zcv = vertex[2]
            vertex_distance = ((xcv - xtp) ** 2 + (ycv - ytp) ** 2) ** 0.5
            # If first checked, automatically gets assigned.
            # If this vertex is nearer than already measured, assign.
            if near_distance is None or vertex_distance < near_distance:
                near_link_id = vertex[0]
                near_vertex_position = vertex[1]
                near_xyz = (xcv, ycv, zcv)
                near_distance = vertex_distance

        # List of link IDs that the routepoint represents should omit the last
        # ID in the walk-back, since that's the last RP's link.
        rp_link_ids = [near_link_id] + candidate_links_tree[near_link_id][:-1]

        self.note('[find_routepoint] Routepoint links %s' % str(rp_link_ids))
        self.note('[find_routepoint] Link vertex %i' % near_vertex_position)

        return rp_link_ids, near_vertex_position, near_xyz


    # Works (2012-05-21).
    def get_od_node_ids(self):
        """Returns the origin/destination node pair IDs as a tuple.

        Requires that the entry & exit nodes be already determined. Currently,
        entry & exit nodes are only added during the fix_routelink_backtracks
        function, so be sure that is run first.

        """

        # No way of determining O&D nodes if there's no links!
        if len(self.routelinks) == 0:
            onode_id, dnode_id = None, None
        else:
            last_position = sorted(self.routelinks)[-1]
            try:
                onode = self.routelinks[1]['entry_node']
                onode_id = onode[0]
            # Means no entry node found on first link.
            except KeyError:
                onode_id = None
            try:
                dnode = self.routelinks[last_position]['exit_node']
                dnode_id = dnode[0]
            # Means no entry node found on final link.
            except KeyError:
                dnode_id = None

        return onode_id, dnode_id


    # Works (2012-05-21).
    def get_routelink_at(self, position):
        """Returns the link ID at the given position in the route."""

        if position in self.routelinks:
            result = self.routelinks[position].get('link_id')
        else:
            result = None

        return result


    # Works (2012-05-25).
    def get_routepoint_at(self, position):
        """Returns a routepoint's info."""

        routepoint_info = {}

        if self.routepoints.get(position):
            for key in self.routepoints[position]:
                routepoint_info[key] = self.routepoints[position].get(key)

        return routepoint_info


    # Works (2012-05-25).
    def get_routepoint_ids(self):
        """Returns all routepoint IDs."""

        routepoint_id_list = []

        for position in self.routepoints:
            routepoint_id_list.append(position)

        return routepoint_id_list


    # Works (2012-05-21).
    def get_routepoints(self):
        """Returns all routepoints along the route."""

        routepoint_list = []

        for position in self.routepoints:
            routepoint_list.append(self.routepoints.get(position))

        return routepoint_list


    # Works (2012-05-21).
    def get_ordered_route_links(self):
        """Returns an ordered list of link IDs."""

        ordered_route_links = []

        for link_position in sorted(set(position for position in self.routelinks)):
            ordered_route_links.append(
                self.routelinks[link_position].get('link_id')
                )

        return ordered_route_links


    # Routes: Data management.

    # Works (2012-07-06): Added multi-link routepoints.
    def add_routepoint_to_route(self, trackpoint_id, link_ids,
                                vertex_position, xyz):
        """Writes link, vertex position & location to given trackpoint.

        This effectively writes to a new trackpoint collection, one which moves
        the points to run along the network.

        """
        self.routepoints[trackpoint_id] = {
            'link_ids': link_ids, 'vertex_position': vertex_position,
            'x': xyz[0], 'y': xyz[1], 'z': xyz[2]
            }

        return trackpoint_id


    # Works (2012-06-29): Added multi-link routepoints.
    def set_entire_route(self, link_network, avoid_problems=True,
                         cursor_size=1):
        """Sets all the routepoints & routelinks using a cursor-like method.

        Cursor size limits the network nodal distance the cursor can check from
        the previous routepoint's location.

        The avoid_problems switch allows the cursor to skip trackpoints that
        may throw the cursor off course. Problem flags are assigned via other
        functions.

        """
        self.note('[set_entire_route] Start.')

        # Initialize the routepoints.
        self.routepoints = {}
        # Add position 1, but with no values in it (needed for later get()).
        self.routepoints[1]= {}

        for trackpoint_id in sorted(self.trackpoints):

            if avoid_problems:
                # Has problems: Skip.
                if self.trackpoints[trackpoint_id].get('problems', False):
                    do_this_one = False
                # Has no problems: Skip.
                else:
                    do_this_one = True
            # Not avoiding, who cares?
            else:
                do_this_one = True

            if do_this_one:
                link_ids, vertex_position, xyz = self.find_routepoint(
                    trackpoint_id, link_network, cursor_size
                    )
                self.add_routepoint_to_route(trackpoint_id, link_ids,
                                             vertex_position, xyz)

        # Align the routelinks with the newly-set routepoints.
        self.set_routelinks()

        # Fix backtrack problems in routelinks.
        ##self.fix_routelink_backtracks(link_network)

        self.note('[set_entire_route] End.')

        return [self.routelinks[link]['link_id'] for link in self.routelinks]


    # Works (2012-07-06): Added multi-link routepoints.
    def set_routelinks(self):
        """Creates the routelink dictionary.

        This can only be effectively run after all routepoints are added to
        the route.

        """

        route_links = []
        link_position = 1
        last_link_id = None

        for routepoint_id in self.routepoints:
            routepoint = self.routepoints.get(routepoint_id)
            # Only add trackpoint info that were calculated into the route.
            if routepoint:
                self.note('[set_routelinks] Routepoint %i.' % routepoint_id)
                # Add link_ids in in travel order (reverse the walk-back).
                for this_link_id in reversed(routepoint['link_ids']):
                    # Only add a link if preceding trackpoint wasn't there also.
                    if this_link_id != last_link_id:
                        route_links.append((link_position, this_link_id))
                        link_position += 1
                        self.note('[set_routelinks] Added link %i.' %
                                  this_link_id)
                    else:
                        self.note('[set_routelinks] Skipping link %i.' %
                                  this_link_id)
                    last_link_id = this_link_id

        # Order the tuples for link positon.
        sorted_route_links = sorted(route_links)

        for link_position, link_id in sorted_route_links:
            self.routelinks[link_position] = {'link_id': link_id}

        return


    # Works (2012-07-09): Added multi-link routepoints; cleaned OD node code.
    def fix_routelink_backtracks(self, link_network):
        """Finds & corrects route links that cause a 'backtrack'.

        If a routepoint is assigned to an incorrect side-link, it will appear
        as if it has 'backtracked' from going down the side-link. In truth,
        the side-link in question was never accessed, it's just a mix-up in the
        assignment. This function will identify these mix-ups, and correct them.

        This function will also add an entry_node & exit_node to routelinks.

        Method of identification: When the route changes links, it cannot enter
        from the same node it entered the previous link with.

        Prescribed fix: Remove side-link from the route links, and re-evaluate
        previous link.

        """
        self.note('[fix_routelink_backtracks] Start.')
        messages = []
        state = True

        # Skip backtracks for no-link: Give O&D nodes NoneType tuples.
        if len(self.routelinks) == 0:
            self.note('[fix_routelink_backtracks] No-link trip.')
            messages.append('No-link trip.')
            self.routelinks[1]['entry_node'] = (None, None)
            self.routelinks[1]['exit_node'] = (None, None)

        # Skip backtracks for single-link: Assign O&D with T&F, respectively.
        elif len(self.routelinks) == 1:
            self.note('[fix_routelink_backtracks] Single-link trip.')
            messages.append('Single-link trip')
            fnode1, tnode1 = link_network.get_link_nodes(self.get_routelink_at(1))
            self.routelinks[1]['entry_node'] = fnode1
            self.routelinks[1]['exit_node'] = tnode1

        # For all other link-counts, correct backtracks & assign true O&Ds.
        else:

            # Start with the first position.
            this_position = 1
            last_position = this_position - 1
            next_position = this_position + 1

            # Loop will run until the current link position doesn't exist.
            while this_position in self.routelinks:

                self.note('[fix_routelink_backtracks] Link position %i.' % (this_position))

                # Get values for current & next link.
                this_link_id = self.get_routelink_at(this_position)
                if last_position in self.routelinks:
                    last_link_id = self.get_routelink_at(last_position)
                else:
                    last_link_id = None
                if next_position in self.routelinks:
                    next_link_id = self.get_routelink_at(next_position)
                else:
                    next_link_id = None

                fnode, tnode = link_network.get_link_nodes(this_link_id)

                # Assign entry & exit nodes.
                # Entry node on from-end.
                if last_link_id in link_network.get_link_ids_at_node(fnode[0], fnode[1]):
                    self.routelinks[this_position]['entry_node'] = fnode
                # Entry node on to-end.
                elif last_link_id in link_network.get_link_ids_at_node(tnode[0], tnode[1]):
                    self.routelinks[this_position]['entry_node'] = tnode
                # last_link_id is None (first position).
                else:
                    self.routelinks[this_position]['entry_node'] = (None, None)
                # Exit node on from-end.
                if next_link_id in link_network.get_link_ids_at_node(fnode[0], fnode[1]):
                    self.routelinks[this_position]['exit_node'] = fnode
                # Exit node on to-end.
                elif next_link_id in link_network.get_link_ids_at_node(tnode[0], tnode[1]):
                    self.routelinks[this_position]['exit_node'] = tnode
                # next_link_id is None (last link position).
                else:
                    self.routelinks[this_position]['exit_node'] = (None, None)

                this_entry_node = self.routelinks[this_position]['entry_node']
                this_exit_node = self.routelinks[this_position]['exit_node']
                if last_position in self.routelinks:
                    last_exit_node = self.routelinks[last_position].get('exit_node')
                else:
                    last_exit_node = (None, None)

                # Check for backtrack.

                # Exit node is the same: This means we've backtracked.
                # Or, last link_id is the same as this one.
                # Ignore if last routelink.
                if (this_exit_node == last_exit_node or
                    this_link_id == last_link_id):
                    # Remove the backtracked link.
                    del self.routelinks[this_position]
                    self.note('[fix_routelink_backtracks] Removed backtrack position %i link %i.' % (this_position, this_link_id))
                    # Reset the routelink order to fill the gap.
                    self.reset_routelink_positions()
                    # Start process over from previous location.
                    this_position -= 1
                # Exit node different: This means a good transition, no backtrack.
                else:
                    # Positive result means we move to the next position.
                    this_position += 1

                last_position = this_position - 1
                next_position = this_position + 1

            # Set origin node.
            try:
                first_position = 1
                first_id = self.get_routelink_at(first_position)
                for node in link_network.get_link_nodes(first_id):
                    if node != self.routelinks[first_position]['exit_node']:
                        self.routelinks[first_position]['entry_node'] = node
            except KeyError:
                self.note('[fix_routelink_backtracks] KeyError while setting o-node.')
                messages.append('KeyError while setting o-node')
                state = False
                try:
                    self.routelinks[first_position]['entry_node'] = (None, None)
                except KeyError:
                    self.note('[fix_routelink_backtracks] KeyError: Backtracks removed all links.')
                    messages.append('KeyError: Backtracks removed all links')

            # Set destination node.
            try:
                final_position = [pos for pos in self.routelinks][-1]
                final_id = self.get_routelink_at(final_position)
                for node in link_network.get_link_nodes(final_id):
                    if node != self.routelinks[final_position]['entry_node']:
                        self.routelinks[final_position]['exit_node'] = node
            except IndexError:
                self.note('[fix_routelink_backtracks] IndexError: Backtracks removed all links.')
                messages.append('IndexError: Backtracks removed all links')
                state = False
            except KeyError:
                self.note('[fix_routelink_backtracks] KeyError while setting d-node.')
                messages.append('KeyError while setting d-node')
                state = False

        self.note('[fix_routelink_backtracks] End.')

        return state, messages


    # Works (2012-05-21).
    def reset_routelink_positions(self):
        """Removes gaps in the link positions (1, 1 + 1, 1 + 2, ...)."""

        old_route = self.routelinks
        new_route = {}
        new_link_position = 1

        for old_link_position in sorted(position for position in old_route):

            new_route[new_link_position] = old_route[old_link_position]
            self.routelinks = new_route
            new_link_position += 1

        return


# Works (2012-07-06): Added link tree with walk-backs.
class LinkNetwork():
    """Stores information for a link network."""

    def __init__(self, source_info={}):
        self.source_path = source_info.get('source_path')
        self.source_data_type = source_info.get('source_data_type')
        self.source_id_field = source_info.get('source_id_field')
        self.source_from_node_field = source_info.get('source_from_node_field')
        self.source_to_node_field = source_info.get('source_to_node_field')
        self.source_from_level_field = source_info.get('source_from_level_field')
        self.source_to_level_field = source_info.get('source_to_level_field')
        self.source_subselect_sql = source_info.get('source_subselect_sql')

        self.link_info = {}
        self.node_info = {}


    # Links: Queries.

    # Testing (2012-07-09): Deprecated max_node_distance.
    def find_near_link_ids(self, link_id):
        """Returns a list of all link IDs near a given link.

        The given link_id is included in the result (that link is near itself).


        """
        nodes_to_cross = list(self.get_link_nodes(link_id))
        near_link_ids = []

        # Get the links at each node.
        for node_id, node_level in nodes_to_cross:
            node_links = self.get_link_ids_at_node(node_id, node_level)
            # Add new link IDs to links_crossed.
            for new_link_id in node_links:
                if new_link_id not in near_link_ids:
                    near_link_ids.append(new_link_id)

        # Remove link ID from input.
        try:
            near_link_ids.remove(link_id)
        except ValueError:
            pass

        return near_link_ids


    # Works (2012-06-29): Deprecate.
    def find_near_link_ids_old(self, link_id, max_node_distance=1):
        """Returns a list of all link IDs near a given link.

        The max_node_distance limits how many nodes can be crossed to reach a
        near-link. If it is 1, only links that share a node with the given link
        will be returned.

        The given link_id is included in the result (that link is near itself).


        """
        crossed_nodes = 0
        # Add values for the 0-link.
        links_crossed = [link_id]
        nodes_to_cross = list(self.get_link_nodes(link_id))
        nodes_crossed = []

        while crossed_nodes < max_node_distance:
            # Get the links at each node.
            for node_id, node_level in nodes_to_cross:
                print('[find_near_link_ids] node = %i.' % node_id)
                node_links = self.get_link_ids_at_node(node_id, node_level)
                print('[find_near_link_ids] node_links = %s.' % str(node_links))

                # Add new link IDs to links_crossed.
                for new_link_id in node_links:
                    print('[find_near_link_ids] new_link_id = %s.' % str(new_link_id))
                    if new_link_id not in links_crossed:
                        print 'not in links_crossed'
                        links_crossed.append(new_link_id)
                    for new_node in self.get_link_nodes(new_link_id):
                        # Add node at other end of new link.
                        if new_node not in nodes_to_cross:
                            nodes_to_cross.append(new_node)
                        else:
                            nodes_to_cross.remove(new_node)

            crossed_nodes += 1

        # Remove link ID from input.
        try:
            links_crossed.remove(link_id)
        except ValueError:
            pass

        return links_crossed


    # Works (unknown date).
    def get_link_count(self):
        """Returns the number of links in the network."""

        return len(link_info)


    # Works (unknown date).
    def get_link_ids(self):
        """Returns a list of all link IDs in the network."""

        return self.link_info.keys()


    # Works (unknown date).
    def get_link_info(self, link_id):
        """Returns all link info for the given link ID."""

        this_link_info = {}
        for info in self.link_info[link_id]:
            this_link_info[info] = self.link_info[link_id].get(info)

        return this_link_info


    # Works (2012-07-06): New code.
    def get_near_link_tree(self, link_id, max_node_distance=1):
        """Returns a tree-dictionary of all link IDs near a given link.

        The tree-dictionary contains all the near link IDs as keys. The values
        for these keys are a list with the walk-back path to the original link.
        (1) The original link will be in the tree, but has an empty walk-back.
        (2) A link adjacent to the original will have only the original link in
        its walk-back.
        (3) Further-out links will have walk-backs larger than one link. All
        links in the walk-back will also be in the tree, and will always end
        with the original link.
        Example: {100: [], 101: [100], 102: [100], 201: [102, 100]}

        The max_node_distance limits how many nodes can be crossed to reach a
        near-link. If it is 1, only links that share a node with the given link
        will be returned.

        If max_node_distance = 0 (no nodes crossed), only the original link
        will be included in the result.

        """
        ##print('[get_near_link_tree] Start.')
        ##print('[get_near_link_tree] Link %i.' % link_id)

        # Add value for the 0-link.
        near_link_tree = {link_id: []}
        ##print('[get_near_link_tree] Added 0-link %i.' % link_id)
        # Set up node collectors.
        crossed_count = 0
        nodes_to_cross = list(self.get_link_nodes(link_id))
        nodes_crossed = []
        ##print('[get_near_link_tree] Nodes to cross %s.' % str(nodes_to_cross))

        while crossed_count < max_node_distance:
            # Evaluate links at each node.
            for node_id, node_level in  sorted(nodes_to_cross):
                link_ids = self.get_link_ids_at_node(node_id, node_level)
                for link_id in link_ids:
                    # Add new links if ID not a key in tree.
                    if link_id not in near_link_tree:
                        ##print('[get_near_link_tree] Link %i not in tree.' % link_id)
                        links_near = self.find_near_link_ids(link_id)
                        ##print link_id, 'links near', links_near
                        for near_link_id in links_near:
                            ##print 'near link', near_link_id
                            # Append the ID & it's tree-chain to the new link.
                            if near_link_id in near_link_tree:
                                near_link_tree[link_id] = (
                                    [near_link_id] + near_link_tree[near_link_id]
                                    )
                                ##print('[get_near_link_tree] Added %i branch %s.' % (link_id, str(near_link_tree[link_id])))
                                break
                    for new_node in self.get_link_nodes(link_id):
                        # Add node at other end of new link.
                        if new_node not in nodes_to_cross and nodes_crossed:
                            nodes_to_cross.append(new_node)
                            ##print('[get_near_link_tree] Added node %s to nodes_to_cross.' % str(new_node))
                nodes_crossed.append((node_id, node_level))
                ##print('[get_near_link_tree] Removed node %i at level %i.' % (node_id, node_level))
            crossed_count += 1
            ##print('[get_near_link_tree] Crossed count = %i.' % crossed_count)

        return near_link_tree


    # Links: Data management.

    # Works (unknown date).
    def add_link_to_network(self, link_id, from_node, to_node,
                            from_level=None, to_level=None, xyz_list=None):
        """Adds a link to the network with the provided information."""

        # Add track ID as key if not there yet.
        self.link_info[link_id] = {
            'from_node': from_node, 'to_node': to_node,
            'from_level': from_level, 'to_level': to_level,
            'points': {}
                }
        position = 1
        if xyz_list:
            for xyz in xyz_list:
                self.link_info[link_id]['points'][position] = {
                    'x': xyz[0], 'y': xyz[1], 'z': xyz[2]
                    }
                position += 1
        self.add_node_to_network(from_node, from_level, [link_id])
        self.add_node_to_network(to_node, to_level, [link_id])

        return


    # Works (unknown date).
    def add_vertex_to_link(self, link_id, xyz, position=None):
        """Adds provided vertex info to a link in the network.

        If position is not noted, will be added at the end.

        """
        if position is None:
            position = len(self.link_info[link_id]['points']) + 2
        self.link_info[link_id]['points'][position] = {
            'x': xyz[0], 'y': xyz[1], 'z': xyz[2]
            }

        return


    # Nodes: Queries.

    # Works (2012-06-29): Added logic to handle incorrect ID or level.
    def get_link_ids_at_node(self, node_id, node_level=None):
        """Returns a list of link IDs for a given node & level.

        If no level is provided, links IDs at all levels of that node will be
        returned.

        """
        # At the least, return an empty list if node_id not there, or
        # node_level not on node.
        link_ids = []

        if node_id in self.node_info:
            if node_level:
                link_ids = self.node_info[node_id].get(node_level)
            else:
                for level in self.node_info[node_id]:
                    link_ids.extend(self.node_info[node_id].get(level))

        return link_ids


    # Works (2012-06-29): Removed list conversion for result (now tuple).
    def get_link_nodes(self, link_id):
        """Returns a pair of tuples containing node IDs & levels."""

        if link_id is None:
            result = ((None, None), (None, None))

        else:
            fnode = (self.link_info[link_id]['from_node'],
                     self.link_info[link_id]['from_level'])
            tnode = (self.link_info[link_id]['to_node'],
                     self.link_info[link_id]['to_level'])
            result = fnode, tnode

        return result


    # Works (2012-06-29).
    def get_node_xyz(self, node_id):
        """Returns a tuple of the XYZ location for a given node."""

        if node_id is None:
            x, y , z = None, None, None
        else:

            node_links = self.get_link_ids_at_node(node_id)

            a_link = node_links[0]

            if node_id == self.link_info[a_link]['from_node']:
                link_position = 1
            elif node_id == self.link_info[a_link]['to_node']:
                link_position = sorted(self.link_info[a_link]['points'])[-1]
            # Recursively impossible, but we'll include it.
            else:
                x, y, z = None, None, None

            x = self.link_info[a_link]['points'][link_position]['x']
            y = self.link_info[a_link]['points'][link_position]['y']
            z = self.link_info[a_link]['points'][link_position]['z']

        return x, y, z


    # Nodes: Data management.

    # Works (unknown date).
    def add_node_to_network(self, node_id, node_level=None, link_ids=[]):
        """Add node ID as key if not there yet. Add level key if not there yet.

        """
        if node_id not in self.node_info:
            self.node_info[node_id] = {}
        if node_level not in self.node_info[node_id]:
            self.node_info[node_id][node_level] = []
        self.node_info[node_id][node_level].extend(link_ids)

        return


    # Vertices: Queries.

    # Works (2012-05-31): New function.
    def get_link_vertex(self, link_id, vertex_position):
        """Returns an xyz tuple for a vertex on a link."""

        # Test for validity of link ID & vertex position.
        if link_id is None:
            link_vertex = (None, None, None)
        else:
            if self.link_info[link_id]['points'].get(vertex_position) is None:
                link_vertex = (None, None, None)

            # All valid, assign tuple data.
            else:
                link_vertex = (
                    self.link_info[link_id]['points'][vertex_position].get('x'),
                    self.link_info[link_id]['points'][vertex_position].get('y'),
                    self.link_info[link_id]['points'][vertex_position].get('z')
                    )

        return link_vertex


    # Works (2012-05-31): Changed to use get_link_vertex().
    def get_link_vertex_tuples(self, link_ids):
        """Returns all the given links' vertices as a list of tuples.

        Format: [(link_id, vertex_position, x, y, z), ...]

        """
        vertex_tuples = []
        for link_id in link_ids:
            for vertex_position in sorted(self.link_info[link_id]['points']):
                vertex_tuples.append(
                    (link_id, vertex_position,
                     self.get_link_vertex(link_id, vertex_position))
                    )

        return vertex_tuples


    # Works (2012-05-31): Changed to use get_link_vertex().
    def get_link_vertices(self, link_id):
        """Returns a list of xyz tuples for all vertices on a link."""

        link_vertices = []

        for vertex_position in sorted(self.link_info[link_id]['points']):
            link_vertices.append(
                self.get_link_vertex(link_id, vertex_position)
                )

        return link_vertices


# Testing (2012-05-07): Not currently used, so unsure of status.
class ModelRoute():
    """Stores information & processes for a model route."""

    def __init__(self, source_info={}):
        self.source_path = source_info.get('source_path')
        self.source_data_type = source_info.get('source_data_type')
        self.source_id_field = source_info.get('source_id_field')
        self.spatial_code = source_info.get('spatial_code')

        self.route_id = source_info.get('route_id')
        self.origin_node = source_info.get('origin_node_id'), None
        self.destination_node = source_info.get('destination_node_id'), None

        self.links = []
        self.nodes = [self.origin_node, self.destination_node]

        # Records the state of the route in processing.
        self.links_sorted = None
        self.nodes_completed = None
        ##self.route_crosses_self = None


    # Route-links: Queries.

    # Works (2012-05-10).
    def get_link_count(self):
        """Returns the number of links in the collection."""

        return len(self.links)


    # Works (2012-05-10).
    def get_link_ids(self):
        """Returns a list of all the link IDs.
         List will be in order only if sort_links_and_nodes() has been done.

         """

        return self.links


    # Route-links: Data management.

    # Testing (2012-05-10).
    def add_link(self, link_id, link_position=None):
        """Inserts a link at the given link_position."""

        # If no link position, just add to the end.
        if link_position is None:
            self.links.append(link_id)
        else:
            self.links.insert(link_position - 1, link_id)

        return self.links


    # Testing (2012-05-10).
    def add_unsorted_links(self, link_ids):
        """Adds a collection of unsorted links."""

        for link_id in link_ids:
            self.add_link(link_id)

        self.links_sorted = None

        return self.links


    # Works (2012-05-10).
    def clear_links(self):
        """Removes all links from the route's inventory."""

        self.links = []

        return self.links


    # Testing (2012-05-10).
    def sort_links_and_nodes(self, link_network):
        """Puts the links in order & adds intermediate nodes."""

        # Reset the nodes list, so we can redo the intermediate nodes.
        self.reset_nodes()

        # Get a list of the link IDs, then clear to write the sort into it.
        unsorted_links = self.get_link_ids()
        self.clear_links()

        # Init markers for iteration & assignment.
        last_link_id = None
        last_node = self.get_node_at(1)

        while unsorted_links:

            possible_link_ids = link_network.get_link_ids_at_node(
                node_id = last_node[0], node_level = last_node[1]
                )
            # Remove the last link ID from the possibles.
            ##if last_link_id is not None:
            ##    possible_link_ids.remove(last_link_id)
            possible_link_ids = [id for id in possible_link_ids if id in unsorted_links]

            # Zero possible means something strange, unknown error.
            if (len(possible_link_ids) == 0 or
                (len(possible_link_ids) == 1 and
                 last_link_id == possible_link_ids[0])):
                self.links_sorted = False
                self.nodes_completed = False
                break

            # One possible, that's the next link.
            # Also: Multiple possibles mean route crosses self. Just pick first.
            else:
                if last_link_id != possible_link_ids[0]:
                    this_link_id = possible_link_ids[0]
                else:
                    this_link_id = possible_link_ids[1]
                self.add_link(this_link_id)
                unsorted_links.remove(this_link_id)

                node01, node02 = link_network.get_link_nodes(this_link_id)

                # Set this_node (as the link's node that's not the last_node).
                if last_node == node01:
                    this_node = node02
                elif last_node == node02:
                    this_node = node01
                # Origin node had no known level: Find node by ID & rewrite o.
                if last_node not in (node01, node02):
                    if last_node[0] == node01[0]:
                        self.reset_nodes(origin_node = node01)
                        last_node = self.origin_node
                        this_node = node02
                    elif last_node[0] == node02[0]:
                        self.reset_nodes(origin_node = node02)
                        last_node = self.origin_node
                        this_node = node01
                self.add_intermediate_node(this_node)
                # Update the link & node markers.
                last_link_id, last_node = this_link_id, this_node

        # After sort iteration, check that last link touches the last two nodes.
        if self.links_sorted is None:
            final_nodes = link_network.get_link_nodes(self.links[-1])
            if (self.get_node_at(1, reverse = True) in final_nodes and
                self.get_node_at(2, reverse = True) in final_nodes):
                self.links_sorted = True
                # Could add a check for node ordering, but not really needed.
                self.nodes_completed = True
                ##self.route_crosses_self = False
            result = self.get_link_ids()

        else:
            result = False

        return result


    # Route-nodes: Queries.

    # Works (2012-05-10).
    def get_node_count(self):
        """Returns the number of nodes in the collection."""

        return len(self.nodes)


    # Testing (unknown date).
    def get_node_at(self, node_position, reverse=False):
        """Returns node tuple at the given position."""

        if reverse:
            result = self.nodes[0 - node_position]
        else:
            result = self.nodes[node_position -1]

        return result


    # Testing (2012-05-10).
    def get_od_nodes(self):
        """Returns the origin/destination node pair IDs in a tuple."""

        return self.get_node_at(1), self.get_node_at(1, reverse = True)


    # Route-nodes: Data management.

    # Testing (2012-05-10).
    def add_intermediate_node(self, node, node_position=None):
        """Inserts a node at the given node position."""

        # If no node position, just add to the end (but before destination).
        if node_position is None:
            self.nodes.insert(-1, node)
        else:
            self.nodes.insert(node_position - 1, node)

        return [node[0] for node in self.nodes]


    # Works (2012-05-10).
    def reset_nodes(self, origin_node=None, destination_node=None):
        """Clear any intermediate nodes, leaving only the O&D nodes.

        If no origin or destination node is provided, will use what's already
        there.

        """
        if origin_node:
            self.origin_node = origin_node
        if destination_node:
            self.destination_node = destination_node

        self.nodes = [self.origin_node, self.destination_node]

        return [node[0] for node in self.nodes]


# Works (2012-06-22).
def create_csv_od_from_arcgis(source_path, data_dictionary):
    """Loads all O&D pairs from an ArcGIS spatial data type to a CSV file.

    Requires a dictionary describing the target data.
    Returns a boolean state and any messages about the process.

    """
    import csv
    try:
        import arcpy
    except ImportError:
        raise Exception

    is_good = True
    messages = []

    # Load O&D nodes from the feature class.
    trip_nodes = {}
    try:
        arc_cursor = arcpy.SearchCursor(source_path)
    except arcpy.ExecuteError:
        is_good = False
        messages.append('Cursor error loading O&D nodes.')

    if is_good:
        for node in arc_cursor:
            trip_id = node.getValue(data_dictionary['trip_id_field'][0])
            node_id = node.getValue(data_dictionary['node_id_field'][0])
            node_type = node.getValue(data_dictionary['node_type_field'][0])
            if trip_id not in trip_nodes:
                trip_nodes[trip_id] = {}
            trip_nodes[trip_id][node_type] = node_id
        try:
            del arc_cursor
        except NameError:
            pass

        # Copy trip_nodes data to CSV output.
        with open(data_dictionary['output_path'], 'wb') as csvfile:
            csv_writer = csv.writer(csvfile)
            # Write the header row.
            csv_writer.writerow([data_dictionary['trip_id_field'][0],
                                 'o' + data_dictionary['node_id_field'][0],
                                 'd' + data_dictionary['node_id_field'][0]])
            # Write the rest of the rows.
            for trip_id in sorted(trip_nodes.keys()):
                csv_writer.writerow([trip_id,
                                     trip_nodes[trip_id]['o'],
                                     trip_nodes[trip_id]['d']])
            # To have a blank row at bottom.
            csv_writer.writerow([])

    return is_good, messages


# Works (2012-06-22).
def create_csv_routes_from_arcgis(source_path, data_dictionary):
    """Loads all route links IDs from an ArcGIS spatial data type to a CSV file.

    Requires a dictionary describing the target data.
    Returns a boolean state and any messages about the process.

    """
    import csv
    try:
        import arcpy
    except ImportError:
        raise Exception

    is_good = True
    messages = []

    # Load route links from the feature class.
    trip_links = {}
    try:
        arc_cursor = arcpy.SearchCursor(source_path)
    except arcpy.ExecuteError:
        is_good = False
        messages.append('Cursor error loading route links.')

    if is_good:
        for link in arc_cursor:
            trip_id = link.getValue(data_dictionary['route_id_field'][0])
            link_id = link.getValue(data_dictionary['link_id_field'][0])
            link_position = link.getValue(
                data_dictionary['link_position_field'][0]
                )
            if trip_id not in trip_links:
                trip_links[trip_id] = {}
            trip_links[trip_id][link_position] = link_id
        try:
            del arc_cursor
        except NameError:
            pass

        # Copy trip_links data to CSV output.
        with open(data_dictionary['output_path'], 'wb') as csvfile:
            csv_writer = csv.writer(csvfile)
            # Write the header row.
            csv_writer.writerow([data_dictionary['route_id_field'][0],
                                 data_dictionary['link_id_field'][0],
                                 data_dictionary['link_position_field'][0]])
            # Write the rest of the rows.
            for trip_id in sorted(trip_links.keys()):
                for link_position in sorted(trip_links[trip_id].keys()):
                    csv_writer.writerow([trip_id,
                                         trip_links[trip_id][link_position],
                                         link_position])
            # To have a blank row at bottom.
            csv_writer.writerow([])

    return is_good, messages


# Works (2012-06-22).
def create_output_layers_for_arcgis(trip_ids, output_dictionary):
    """Creates an ArcGIS group layer for each of the trip IDs provided.

    This function will replace the group layers (if they exist) of the provided
    trips with a newer one.

    Requires one or more trip IDs and a dictionary describing the target data.

    """
    import os
    try:
        import arcpy
        from arcpy.mapping import Layer
    except ImportError:
        raise Exception

    is_good = True
    folder_path = output_dictionary['folder_path']
    name_format = output_dictionary['name_format']
    template_path = output_dictionary['template_path']
    sublayers_sql = output_dictionary['sublayers_sql']
    messages = []

    while is_good:

        # Load the template_layer
        template_layer = Layer(template_path)

        for trip_id in trip_ids:

            trip_layer_name = name_format.replace('<ID>', str(trip_id))
            trip_layer_path = os.path.join(
                folder_path, '%s.lyr' % trip_layer_name
                )

            # Remove old layer for trip.
            try:
                if arcpy.Exists(trip_layer_path):
                    arcpy.management.Delete(trip_layer_path)
            except arcpy.ExecuteError:
                is_good = False
                messages.append(
                    'Error: Create output layers failed in removing old layer.'
                    )
                break

            # Apply the new layer name.
            template_layer.name = trip_layer_name

            # Apply the new definition queries.
            for sublayer in template_layer:
                sublayer.definitionQuery = sublayers_sql.get(
                    sublayer.name, '').replace('<ID>', str(trip_id))

            # Save a copy under the new name.
            template_layer.saveACopy(trip_layer_path)

        del template_layer
        break

    return is_good, messages


# Works, takes forever (2012-06-22).
def get_all_trip_ids(trip_path, trip_id_field):
    """Returns a list of all trip IDs that exist in a given trip collection."""

    try:
        import arcpy
    except ImportError:
        raise Exception

    trip_ids = []
    trip_cursor = arcpy.SearchCursor(trip_path)
    for trip_part in trip_cursor:
        if trip_part.getValue(trip_id_field) is not None:
            if trip_part.getValue(trip_id_field) not in trip_ids:
                trip_ids.append(trip_part.getValue(trip_id_field))

    return trip_ids


# Works (2012-05-21).
def load_gps_trips_from_mysql(data_dictionary, trip_ids=None,
                              logfile_path=None, trip_limit=None):
    """Loads a given track from an MySQL database.

    Requires the MySQLdb library to be installed, and a conversion
    dictionary for the source data. Returns a list of the GPSTrip class
    instances for the given trip IDs (optional: omitting gets all in DB) and
    omitting all trips noted in the logfile (optional).

    """
    import csv
    try:
        import MySQLdb
    except ImportError:
        raise Exception

    table_name = data_dictionary.get('source_table')
    trip_id_field = data_dictionary.get('source_id_field')
    timestamp_field = data_dictionary.get('source_timestamp_field')
    easting_field = data_dictionary.get('source_easting_field')
    northing_field = data_dictionary.get('source_northing_field')
    elevation_field = data_dictionary.get('source_elevation_field')
    xy_accuracy_field = data_dictionary.get('source_xy_accuracy_field')
    z_accuracy_field = data_dictionary.get('source_z_accuracy_field')

    # Connect to DB.
    conn = MySQLdb.Connect(host = data_dictionary.get('source_host'),
                           port = data_dictionary.get('source_port'),
                           user = data_dictionary.get('source_user'),
                           passwd = data_dictionary.get('source_password'),
                           db = data_dictionary.get('source_db'),
                           compress = 1)

    # Gather all wanted trip IDs from the DB.
    trips_from_database = []
    cursor = conn.cursor(cursorclass = MySQLdb.cursors.DictCursor)
    cursor.execute("select distinct %s from %s" % (trip_id_field, table_name))
    gps_trip_ids = cursor.fetchall()
    for row in gps_trip_ids:
        # Including all trips.
        if trip_ids is None:
            trips_from_database.append(int(row[trip_id_field]))
        # Including just the given trips.
        else:
            if int(row[trip_id_field]) in trip_ids:
                trips_from_database.append(int(row[trip_id_field]))

    trips_to_analyze = []

    # Check the logfile (if included).
    if logfile_path:
        # Gather all logged trip IDs.
        trips_logged = []
        log_rows = csv.DictReader(open(logfile_path, 'rU'))
        for row in log_rows:
            trip_logged_as = row.get('analysis_complete', str(None))
            if trip_logged_as in (str(True), str(False)):
                trips_logged.append(int(row.get(trip_id_field)))
        # Omit logged trips.
        for trip_id in trips_from_database:
            if trip_id in trips_logged:
                pass
            else:
                trips_to_analyze.append(trip_id)
    else:
        trips_to_analyze = trips_from_database

    # If a trip limit is set, truncate trips_to_analyze.
    if trip_limit:
        trips_to_analyze = trips_to_analyze[:trip_limit]

    # Gather all rows for the chosen trips.
    all_trips_rows = []
    for trip_id in  trips_to_analyze:
        cursor.execute(
            "select * from %s where trip_id = %s order by %s, %s" %
            (table_name, str(trip_id), trip_id_field, timestamp_field)
            )
        trip_rows = cursor.fetchall()
        all_trips_rows.append(trip_rows)

    # Close MySQL connection.
    cursor.close()
    conn.close()

    # Create class instances and write trackpoints into it.
    all_trips = []
    for trip_rows in all_trips_rows:
        if trip_rows:
            data_dictionary['trip_id'] = trip_rows[0].get(trip_id_field)
            this_trip = GPSTrip(data_dictionary)
            track_position = 1
            # Add each trackpoint to the trip class instance.
            for row in trip_rows:
                # Convert GPS lat/lon to state plane feet.
                x, y = transform_gps_to_orsps(
                    longitude = row.get(easting_field),
                    latitude = row.get(northing_field)
                    )
                this_trip.add_trackpoint(
                    position = track_position,
                    timestamp = row.get(timestamp_field),
                    xyz = (x, y, row.get(elevation_field)),
                    accuracy = (row.get(xy_accuracy_field),
                                row.get(z_accuracy_field))
                    )
                track_position += 1

            all_trips.append(this_trip)

    return all_trips


# Works (2012-05-07).
def load_link_network_from_arcgis(data_dictionary):
    """Loads the link network from an ArcGIS-native spatial data type.

    Requires a conversion dictionary for the source data. Returns a dictionary
    containing each individual link as a sub-dictionary.

    """
    try:
        import arcpy
    except ImportError:
        raise Exception

    network = LinkNetwork(data_dictionary)

    try:
        ncursor = arcpy.SearchCursor(network.source_path,
                                     network.source_subselect_sql)
    except arcpy.ExecuteError:
        print 'Need to add a error prompt here.'

    for n in ncursor:
        if network.source_from_level_field:
            fzlev = n.getValue(network.source_from_level_field)
        if network.source_to_level_field:
            tzlev = n.getValue(network.source_to_level_field)

        network.add_link_to_network(
            link_id = n.getValue(network.source_id_field),
            from_node = n.getValue(network.source_from_node_field),
            to_node = n.getValue(network.source_to_node_field),
            from_level = fzlev,
            to_level = tzlev,
            xyz_list = None
            )

        # Get the points from the link's geometry (only set for single-part).
        point_position = 1
        for point in n.getValue(data_dictionary['shape_field']).getPart(0):
            network.add_vertex_to_link(
                link_id = n.getValue(network.source_id_field),
                xyz = (point.X, point.Y, point.Z), position = point_position
                )
            point_position += 1

    del ncursor

    return network


# Works, unused (2012-05-21).
def load_model_routes_from_csv_folder(data_dictionary):
    """Loads model routes from a folder full of CSVs."""

    import csv
    import os

    # Get list of files in CSV folder & filter for correct filenames.
    folder = data_dictionary['source_path']
    files_in_folder = os.listdir(folder)
    files_to_load = []
    for filename in files_in_folder:
        if filename.endswith('.csv'):
            files_to_load.append(filename)

    model_routes = []

    # Parse name & create a model route instance.
    for filename in files_to_load:

        # For now, we're going to assume the format & order in the filename.
        # We'll get sneaky with regular expressions later.
        import re
        id_parts = re.split(r'[\D]', filename)
        while '' in id_parts:
            id_parts.remove('')
        route_id, origin_node_id, destination_node_id = [int(part) for part in id_parts]

        this_route_data = data_dictionary
        this_route_data['source_path'] = os.path.join(folder, filename)
        this_route_data['route_id'] = route_id
        this_route_data['origin_node_id'] = origin_node_id
        this_route_data['destination_node_id'] = destination_node_id

        # Load the route links from the CSV file.
        csv_rows = csv.DictReader(open(this_route_data['source_path'], 'rU'))
        unsorted_links = [int(row[data_dictionary['source_link_id_field']]) for row in csv_rows]

        # Check whether CSV has any link-rows before bothering with it.
        if unsorted_links:

            # Construct this_route as a class.
            this_route = ModelRoute(this_route_data)
            this_route.add_unsorted_links(unsorted_links)

            # Add this_route to the list of routes.
            model_routes.append(this_route)

    return model_routes


# Works (2012-05-21).
def transform_gps_to_orsps(longitude, latitude):
    """Converts GPS latitude & longitude into our in-house coordinate system.

    Specifically, this converts them into Oregon State Plane South, NAD 83 HARN
    in international feet. Returns a tuple of the new x/y coordinates.

    This function requires the local computer have the pyproj library installed.

    """
    try:
        from pyproj import Proj
    except ImportError:
        raise Exception

    sps_proj = Proj(proj = 'lcc',  # Lambert conformal conic.
                    lat_1 = 42.33333333333334, lat_2 = 44,
                    lat_0 = 41.66666666666666, lon_0 = -120.5,
                    x_0 = 1500000, y_0 = 0,
                    ellps = 'GRS80', units = 'ft',
                    no_defs = True, preserve_units = True)

    x, y = sps_proj(longitude, latitude)

    return x, y


# Works (2012-05-25).
def update_od_nodes_in_arcgis(gps_trips, link_network, output_dictionary):
    """Update an ArcGIS-native spatial class of the origin & destination points.

    This function will replace old features (if they exist) of the provided
    GPS trips in the output data with what is current in the GPSTrips class.

    Requires one or more GPSTrip classes, the link network class, and a
    dictionary describing the target data.

    """
    try:
        import arcpy
        from arcpy import management as mgmt
    except ImportError:
        raise Exception

    oob_points = []

    is_good = True
    output_path = output_dictionary['output_path']
    event_id_field = output_dictionary['trip_id_field'][0]
    node_id_field = output_dictionary['node_id_field'][0]
    node_type_field = output_dictionary['node_type_field'][0]
    shape_field = output_dictionary['shape_field']
    messages = []

    while is_good:

        # Collect all trip_ids.
        trip_ids = []
        for trip in gps_trips:
            trip_ids.append(trip.trip_id)

        # Remove features for all the trips in question.
        try:
            layer_sql = "%s in (%s)" % (
                event_id_field, ', '.join([str(tid) for tid in trip_ids])
                )
            mgmt.MakeFeatureLayer(in_features = output_path,
                                  out_layer = 'output_layer',
                                  where_clause = layer_sql)
            mgmt.DeleteFeatures('output_layer')
            mgmt.Delete('output_layer')
        except arcpy.ExecuteError:
            is_good = False
            messages.append(
                'Error: Update O&D pairs failed in removing old trip rows.'
                )
            break

        try:
            output_cursor = arcpy.InsertCursor(output_path)
        except arcpy.ExecuteError:
            is_good = False
            messages.append(
                'Error: Update O&D pairs failed to create a cursor.'
                )
            break

        # Iterate through all the trips.
        for trip in gps_trips:
            origin_id, destination_id = trip.get_od_node_ids()

            for node_type in ('o', 'd'):
                if node_type == 'o':
                    node_id = origin_id
                else:
                    node_id = destination_id
                node_xyz = link_network.get_node_xyz(node_id)

                # Make new feature & assign attributes/geometry.
                out_row = output_cursor.newRow()

                # Build point-shape from node_xyz.
                # (X, Y, Z, M, ID)
                point = arcpy.Point(X = node_xyz[0],
                                    Y = node_xyz[1],
                                    Z = node_xyz[2])
                try:
                    out_row.setValue(shape_field, point)
                # Out-of-bounds geometry throws a RuntimeError.
                except RuntimeError:
                    state = False
                    oob_points.append('trip %s-node %s' %
                                      (str(trip.trip_id), str(node_id)))

                out_row.setValue(event_id_field, trip.trip_id)
                out_row.setValue(node_id_field, node_id)
                out_row.setValue(node_type_field, node_type)

                output_cursor.insertRow(out_row)

        try:
            del output_cursor
            del out_row
        except NameError:
            pass
        break

    # Build out-of-bounds message.
    if oob_points:
        messages.append('Out-of-bounds node geometry (%s)' %
                        '; '.join(oob_points))

    return is_good, messages


# Works (2012-05-25).
def update_routes_in_arcgis(gps_trips, link_network, output_dictionary):
    """Updates an ArcGIS-native spatial class  of routes (linkwise).

    This function will replace old features (if they exist) of the provided
    GPS trips in the output data with what is current in the GPSTrips class.

    Requires one or more GPSTrip classes, a link network, and a dictionary
    describing the target data.

    """
    try:
        import arcpy
        from arcpy import management as mgmt
    except ImportError:
        raise Exception

    is_good = True
    output_path = output_dictionary['output_path']
    event_id_field = output_dictionary['route_id_field'][0]
    link_id_field = output_dictionary['link_id_field'][0]
    link_position_field = output_dictionary['link_position_field'][0]
    shape_field = output_dictionary['shape_field']
    messages = []

    while is_good:

        # Collect all trip_ids.
        trip_ids = []
        for trip in gps_trips:
            trip_ids.append(trip.trip_id)

        # Remove features for all the trips being updated.
        try:
            layer_sql = "%s in (%s)" % (
                event_id_field, ', '.join([str(tid) for tid in trip_ids])
                )
            mgmt.MakeFeatureLayer(in_features = output_path,
                                  out_layer = 'output_layer',
                                  where_clause = layer_sql)
            mgmt.DeleteFeatures('output_layer')
            mgmt.Delete('output_layer')
        except arcpy.ExecuteError:
            is_good = False
            messages.append(
                'Error: Update routes failed in removing old trip rows.'
                )
            break

        try:
            output_cursor = arcpy.InsertCursor(output_path)
        except arcpy.ExecuteError:
            is_good = False
            messages.append('Error: Update routes failed to create a cursor.')
            break

        # Iterate through all the trips.
        for trip in gps_trips:

            link_position = 1

            for link_id in trip.get_ordered_route_links():
                # Make new feature & assign attributes/geometry.
                out_row = output_cursor.newRow()
                link_vertices = link_network.get_link_vertices(link_id)
                # Build line-shape from routepoints.
                line_array = arcpy.Array()

                for vertex in link_vertices:
                    # (X, Y, Z, M, ID)
                    point = arcpy.Point(X = vertex[0],
                                        Y = vertex[1],
                                        Z = vertex[2])
                    line_array.add(point)

                out_row.setValue(shape_field, line_array)
                out_row.setValue(event_id_field, trip.trip_id)
                out_row.setValue(link_id_field, link_id)
                out_row.setValue(link_position_field, link_position)

                output_cursor.insertRow(out_row)
                link_position += 1


        try:
            del output_cursor
            del out_row
        except NameError:
            pass
        break

    return is_good, messages


# Works (2012-05-25).
def update_routepoints_in_arcgis(gps_trips, output_dictionary):
    """Update an ArcGIS-native spatial class of routepoints.

    This function will replace old features (if they exist) of the provided
    GPS trips in the output data with what is current in the GPSTrips class.

    Requires one or more GPSTrip classes and a dictionary describing the target
    data.

    """
    try:
        import arcpy
        from arcpy import management as mgmt
    except ImportError:
        raise Exception

    oob_points = []

    is_good = True
    output_path = output_dictionary['output_path']
    event_id_field = output_dictionary['route_id_field'][0]
    routepoint_id_field = output_dictionary['routepoint_id_field'][0]
    link_ids_field = output_dictionary.get('link_ids_field', (None,))[0]
    link_vertex_field = output_dictionary.get('link_vertex_field', (None,))[0]
    shape_field = output_dictionary['shape_field']
    messages = []

    while is_good:

        # Collect all trip_ids.
        trip_ids = []
        for trip in gps_trips:
            trip_ids.append(trip.trip_id)

        # Remove features for all the trips being updated.
        try:
            layer_sql = "%s in (%s)" % (
                event_id_field, ', '.join([str(tid) for tid in trip_ids])
                )
            mgmt.MakeFeatureLayer(in_features = output_path,
                                  out_layer = 'output_layer',
                                  where_clause = layer_sql)
            mgmt.DeleteFeatures('output_layer')
            mgmt.Delete('output_layer')
        except arcpy.ExecuteError:
            is_good = False
            messages.append(
                'Error: Update routepoints failed in removing old trip rows.'
                )
            break
        try:
            output_cursor = arcpy.InsertCursor(output_path)
        except arcpy.ExecuteError:
            is_good = False
            messages.append(
                'Error: Update routepoints failed to create a cursor.'
                )
            break

        # Iterate through all the trips.
        for trip in gps_trips:

            for routepoint_id in trip.get_routepoint_ids():
                routepoint = trip.get_routepoint_at(routepoint_id)

                # Make new feature & assign attributes/geometry.
                out_row = output_cursor.newRow()

                # Build point-shape from routepoint.
                # (X, Y, Z, M, ID)
                point = arcpy.Point(X = routepoint['x'],
                                    Y = routepoint['y'],
                                    Z = routepoint['z'])

                try:
                    out_row.setValue(shape_field, point)
                # Out-of-bounds geometry throws a RuntimeError.
                except RuntimeError:
                    oob_points.append('trip %s-routepoint %s' %
                                      (str(trip.trip_id), str(routepoint_id)))

                out_row.setValue(event_id_field, trip.trip_id)
                out_row.setValue(routepoint_id_field, routepoint_id)

                if link_ids_field:
                    link_ids = routepoint.get('link_ids', [])
                    link_ids_string = ''
                    for link_id in link_ids:
                        if len(link_ids_string) == 0:
                            new_string = str(link_id)
                        else:
                            new_string = ','.join([link_ids_string, str(link_id)])
                        if len(new_string) <= 64:
                            link_ids_string = new_string
                        else:
                            break
                    out_row.setValue(link_ids_field, link_ids_string)

                if link_vertex_field:
                    out_row.setValue(link_vertex_field,
                                     routepoint.get('vertex_position'))

                output_cursor.insertRow(out_row)

        try:
            del output_cursor
            del out_row
        except NameError:
            pass
        break

    # Build out-of-bounds message.
    if oob_points:
        messages.append('Out-of-bounds routepoint geometry (%s)' %
                        '; '.join(oob_points))

    return is_good, messages


# Works (2012-05-25).
def update_trackpoints_in_arcgis(gps_trips, output_dictionary):
    """Update an ArcGIS-native spatial class of trackpoints.

    This function will replace old features (if they exist) of the provided
    GPS trips in the output data with what is current in the GPSTrips class.

    Requires one or more GPSTrip classes and a dictionary describing the target
    data.

    """
    try:
        import arcpy
        from arcpy import management as mgmt
    except ImportError:
        raise Exception

    oob_points = []

    is_good = True
    output_path = output_dictionary['output_path']
    event_id_field = output_dictionary['track_id_field'][0]
    trackpoint_id_field = output_dictionary['trackpoint_id_field'][0]
    link_ids_field = output_dictionary.get('link_ids_field', (None, None))[0]
    link_vertex_field = output_dictionary.get('link_vertex_field',
                                              (None, None))[0]
    shape_field = output_dictionary['shape_field']
    messages = []

    while is_good:

        # Collect all trip_ids.
        trip_ids = []
        for trip in gps_trips:
            trip_ids.append(trip.trip_id)

        # Remove features for all the trips being updated.
        try:
            layer_sql = "%s in (%s)" % (
                event_id_field, ', '.join([str(tid) for tid in trip_ids])
                )
            mgmt.MakeFeatureLayer(in_features = output_path,
                                  out_layer = 'output_layer',
                                  where_clause = layer_sql)
            mgmt.DeleteFeatures('output_layer')
            mgmt.Delete('output_layer')
        except arcpy.ExecuteError:
            is_good = False
            messages.append(
                'Error: Update trackpoints failed in removing old trip rows.'
                )
            break

        try:
            output_cursor = arcpy.InsertCursor(output_path)
        except arcpy.ExecuteError:
            is_good = False
            messages.append(
                'Error: Update trackpoints failed to create a cursor.'
                )
            break

        # Iterate through all the trips.
        for trip in gps_trips:

            for trackpoint_id in trip.get_trackpoint_ids():
                trackpoint = trip.get_trackpoint_at(trackpoint_id)
                # Make new feature & assign attributes/geometry.
                out_row = output_cursor.newRow()
                # Build point-shape from trackpoint.
                # (X, Y, Z, M, ID)
                point = arcpy.Point(X = trackpoint['x'],
                                    Y = trackpoint['y'],
                                    Z = trackpoint['z'])
                try:
                    out_row.setValue(shape_field, point)
                # Out-of-bounds geometry throws a RuntimeError.
                except RuntimeError:
                    state = False
                    oob_points.append('trip %s-trackpoint %s' %
                                      (str(trip.trip_id), str(trackpoint_id)))

                out_row.setValue(event_id_field, trip.trip_id)
                out_row.setValue(trackpoint_id_field, trackpoint_id)

                routepoint = trip.get_routepoint_at(trackpoint_id)
                if link_ids_field:
                    link_id_strings = [str(lid) for lid in
                                       routepoint.get('link_ids', [])]
                    links_string = ';'.join(link_id_strings)[:64]
                    out_row.setValue(link_ids_field, links_string)
                if link_vertex_field:
                    out_row.setValue(link_vertex_field,
                                     routepoint.get('vertex_position'))

                output_cursor.insertRow(out_row)

        try:
            del output_cursor
            del out_row
        except NameError:
            pass
        break

    # Build out-of-bounds message.
    if oob_points:
        messages.append('Out-of-bounds trackpoint geometry (%s)' %
                        '; '.join(oob_points))

    return is_good, messages
