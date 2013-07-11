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
##DEMO_GPS_TRIP_FROM_MYSQL = {
##    'source_host': 'database.domain-name-string.org',
##    'source_port': 3306,  # Port as integer.
##    'source_user': 'username_string',
##    'source_password': 'user_password_string',
##    'source_db': 'database_name_string',
##    'source_table': 'table_name_as_string',
##    'source_data_type': 'MySQL',  # Not yet used; just for info.
##    'source_id_field': 'id_field_string',
##    'source_timestamp_field': 'timestamp_field_name_string',
##    'source_easting_field': 'easting_field_name_string',
##    'source_northing_field': 'northing_field_name_string',
##    'source_elevation_field': 'elevation_field_name_string',
##    'source_xy_accuracy_field': 'xy_accuracy_field_name_string',
##    'source_z_accuracy_field': 'z_accuracy_field_name_string',
##    }
##
##
### Works (2012-07-06): Added multi-link routepoints; added debug notes.
##class GPSTrip():
##    """Stores trackpoint & other information for a GPS-recorded trip."""
##
##
##    # Works (2012-05-21).
##    def __init__(self, source_info={}):
##
##        self.trip_id = source_info.get('trip_id')
##
##        self.source_path = source_info.get('source_path')
##        self.source_data_type = source_info.get('source_data_type')
##        self.source_id_field = source_info.get('source_id_field')
##        self.source_timestamp_field = source_info.get('source_timestamp_field')
##        self.source_elevation_field = source_info.get('source_elevation_field')
##        self.source_xy_accuracy_field = source_info.get('source_xy_accuracy_field')
##        self.source_z_accuracy_field = source_info.get('source_z_accuracy_field')
##
##        self.statistics = {}
##
##        self.trackpoints = {}
##        self.track_speeds = {}
##        self.routepoints = {}
##        self.routelinks = {}
##
##        self.all_notes = []  # Repository for debugging notes.
##        self.echo = False  # Will print debugging notes to the console.
##        self.note('[__init__] trip_id %i.' % self.trip_id)
##
##    def note(self, note_to_add):
##        """Adds a note to self.all_notes."""
##
##        timestamp = str(datetime.datetime.now().replace(microsecond = 0))
##
##        self.all_notes.append('%s: %s' % (timestamp, note_to_add))
##
##        if self.echo:
##            print note_to_add
##
##        return
##
##
##    # Tracks: Queries.
##
##    # Works (2012-07-09).
##    def find_absolute_trip_length(self):
##        """Calculates the absolute length of the trip.
##
##        Absolute distance is a pure calculation of distance traveled between
##        the ordered trackpoints in the trip. This is *not* the distance
##        along the route as assigned by the map-matching.
##
##        Returns a tuple containing:
##        (absolute_distance, absolute_time_seconds)
##
##        """
##
##        # Init travel tracker.
##        distance_traveled = 0.0
##        time_traveled_seconds = 0.0
##
##        for this_id in self.trackpoints:
##
##            if this_id != 1:
##                that_id = this_id - 1
##
##                # This trackpoint's geometry.
##                this_x = self.trackpoints[this_id].get('x')
##                this_y = self.trackpoints[this_id].get('y')
##                this_time = self.trackpoints[this_id].get('timestamp')
##
##                # Previous trackpoint's geometry.
##                that_x = self.trackpoints[that_id].get('x')
##                that_y = self.trackpoints[that_id].get('y')
##                that_time = self.trackpoints[that_id].get('timestamp')
##
##                # Find distance & time between points.
##                this_distance = (
##                    ((that_x - this_x) ** 2 + (that_y - this_y) ** 2) ** 0.5
##                    )
##                this_time_seconds = (this_time - that_time).seconds
##
##                # Add distance & time to total.
##                distance_traveled += this_distance
##                time_traveled_seconds += this_time_seconds
##
##        return distance_traveled, time_traveled_seconds
##
##
##    # Works (2012-05-31).
##    def find_largest_trackpoint_gap(self):
##        """Determines the largest gap between adjacent trackpoints.
##
##        Returns a tuple containing:
##        (largest_gap_distance, largest_gap_time_seconds)
##
##        """
##
##        # Init gap tracker.
##        largest_gap_distance = 0.0
##        largest_gap_time_seconds = 0.0
##
##        for this_id in self.trackpoints:
##
##            if this_id != 1:
##                that_id = this_id - 1
##
##                # This trackpoint's geometry.
##                this_x = self.trackpoints[this_id].get('x')
##                this_y = self.trackpoints[this_id].get('y')
##                this_time = self.trackpoints[this_id].get('timestamp')
##
##                # Previous trackpoint's geometry.
##                that_x = self.trackpoints[that_id].get('x')
##                that_y = self.trackpoints[that_id].get('y')
##                that_time = self.trackpoints[that_id].get('timestamp')
##
##                this_gap_distance = (
##                    ((that_x - this_x) ** 2 + (that_y - this_y) ** 2) ** 0.5
##                    )
##
##                # If gap_distance largest sor far, assign.
##                if this_gap_distance > largest_gap_distance:
##                    largest_gap_distance = this_gap_distance
##                    largest_gap_time_seconds = (this_time - that_time).seconds
##
##        return largest_gap_distance, largest_gap_time_seconds
##
##
##    # Doesn't currently work (2012-07-09).
##    def find_trackpoints_outside_area(self, easting_range, northing_range):
##        """Compiles all trackpoints that fall outside the given area.
##
##        The easting & northing ranges represent an envelope that defines the
##        "area".
##
##        Returns a list of trackpoint IDs that are outside the envelope.
##
##        """
##
##        # Init ID list.
##        trackpoints_outside = []
##
##        for trackpoint_id in self.trackpoints:
##
##            # This trackpoint's geometry.
##            x = self.trackpoints[trackpoint_id].get('x')
##            y = self.trackpoints[trackpoint_id].get('y')
##
##            x_is_outside = x < easting_range[0] or x > easting_range[1]
##            y_is_outside = y < northing_range[0] or y > northing_range[1]
##
##            if x_is_outside or y_is_outside:
##                trackpoints_outside.append(trackpoint_id)
##
##        return trackpoints_outside
##
##
##    # Works (2012-05-21).
##    def get_trackpoint_count(self):
##        """Returns the number of trackpoint in the trip."""
##
##        return len(self.trackpoints)
##
##
##    # Works (2012-05-25).
##    def get_trackpoint_ids(self):
##        """Returns all trackpoint IDs."""
##
##        trackpoint_id_list = []
##
##        for position in self.trackpoints:
##            trackpoint_id_list.append(position)
##
##        return trackpoint_id_list
##
##
##    # Works (2012-05-25).
##    def get_trackpoint_at(self, position):
##        """Returns a trackpoint's info."""
##
##        trackpoint_info = {}
##
##        for key in self.trackpoints.get(position):
##            trackpoint_info[key] = self.trackpoints[position].get(key)
##
##        return trackpoint_info
##
##
##    # Works (2012-05-25).
##    def get_trackpoints(self):
##        """Returns all trackpoint info."""
##
##        trackpoint_list = []
##
##        for position in self.trackpoints:
##            trackpoint_list.append(self.trackpoints.get(position))
##
##        return trackpoint_list
##
##
##    # Tracks: Data management.
##
##    # Works (2012-05-21).
##    def add_trackpoint(self, position, timestamp, xyz, accuracy):
##        """Adds a single trackpoint to the trip's dictionary"""
##
##        self.trackpoints[position] = {
##            'trackpoint_id': position, 'timestamp': timestamp,
##            'x': xyz[0], 'y': xyz[1], 'z': xyz[2],
##            'xy_accuracy': accuracy[0], 'z_accuracy': accuracy[1]
##            }
##
##        return position, xyz
##
##
##    # Not sure this works at all (2012-05-21).
##    def identify_problem_trackpoints(self):
##        """Identifies & marks suspicious trackpoints via specific QA methods."""
##
##        speed_mean = self.statistics['speed_mean']
##        speed_stdev = self.statistics['speed_standard_deviation']
##
##        for position in sorted(self.trackpoints):
##            # Initial position is considered gospel. If bad, need to ID before
##            # running most of these class functions (and map-matching).
##            if position > 1:
##                ##lastposition = position - 1
##                thispoint = self.trackpoints[position]
##                ##lastpoint = self.trackpoints[lastposition]
##                # Move another point back to until no previous problems, or
##                # you've reached the initial point (which is gospel).
##                ##while lastpoint.get('problems') and lastposition != 1:
##                ##    lastposition = lastposition - 1
##                ##    lastpoint = self.trackpoints[lastposition]
##
##                # Problem: Unlikely speed.
##                # If speed is more than two standard deviations from the mean,
##                # mark.
##                deviation = abs(speed_mean - self.track_speeds['speed'])
##                if deviation >= speed_stdev * 2:
##                    if thispoint['problems']:
##                        thispoint['problems'].append('speed')
##                    else:
##                        thispoint['problems'] = ['speed']
##
##        return
##
##
##    # Not yet written (2012-05-21).
##    def reset_trackpoint_positions(self):
##        """Resets trackpoint position values based on timestamp values.
##
##        Still needs to be written.
##
##        """
##
##        return
##
##
##    # Not working, unused (2012-005-21).
##    def set_trackpoint_speed(self):
##        """Assigns a speed to a trackpoint based on the previous one."""
##
##        for position in self.trackpoints:
##         # Initial position gets a speed of 0 (no previous point, plus logic
##         # says you aren't moving until you've begun).
##         if position == 1:
##            self.track_speeds[position] = 0.0
##         else:
##            this_time = self.trackpoints[position]['timestamp']
##            that_time = self.trackpoints[position - 1]['timestamp']
##            dt = this_time - that_time
##            time_elapsed = dt.seconds + dt.microseconds * 10 ** -6
##            x2 = self.trackpoints[position]['x']
##            y2 = self.trackpoints[position]['y']
##            x1 = self.trackpoints[position - 1]['x']
##            y1 = self.trackpoints[position - 1]['y']
##            distance_elapsed = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
##
##            self.track_speeds[position] = (distance_elapsed / time_elapsed)
##
##        # Update speed statistics.
##
##        # Mean speed.
##        speed_sum = 0
##        speed_count = 0
##        for position in self.trackpoints:
##            # Skip the initial point for speed stats; logic precludes this.
##            if position != 1:
##                speed_sum += self.track_speeds[position]
##                speed_count += 1
##        speed_mean = speed_sum / speed_count
##
##        # Standard deviation.
##        for position in self.trackpoints:
##            # Skip the initial point for speed stats; logic precludes this.
##            variance_numer = 0
##            variance_denom = 0
##            if position != 1:
##                variance_numer += (
##                    self.track_speeds[position] - speed_mean) ** 2
##                variance_denom += 1
##        speed_variance = variance_numer / variance_denom
##        speed_stdev = variance ** 0.5
##
##        self.statistics['speed_mean'] = speed_mean
##        self.statistics['speed_standard_deviation'] = speed_stdev
##
##        return speed_mean
##
##
##    # Routes: Queries.
##
##    # Works (2012-07-06): Added multi-link routepoints.
##    def find_routepoint(self, trackpoint_id, link_network, cursor_size=1):
##        """Finds the nearest vertex in the link network.
##
##        Previous routepoint must have been already solved.
##
##        """
##        self.note('[find_routepoint] Trackpoint %i.' % trackpoint_id)
##        self.note('[find_routepoint] Info %s.' % str(self.get_trackpoint_at(trackpoint_id)))
##
##        # Update potential link IDs.
##
##        # If initial trackpoint, all links are fair game. Same goes if
##        # something else weird happens (in that case, link_ids will be [None]).
##        if trackpoint_id == 1:
##            self.note('[find_routepoint] Is first trackpoint.')
##            candidate_links = link_network.get_link_ids()
##            candidate_links_tree = {}
##            for candidate_link_id in candidate_links:
##                candidate_links_tree[candidate_link_id] = []
##        # Otherwise, get only links sharing a node with the cursor link.
##        else:
##            last_routepoint_id = trackpoint_id - 1
##            last_link_ids = self.routepoints[last_routepoint_id].get('link_ids', [])
##            while len(last_link_ids) == 0:
##                # You're back at the startpoint and have no link_id. All links.
##                if last_routepoint_id == 0:
##                    self.note('[find_routepoint] Back to first routepoint.')
##                    candidate_links = link_network.get_link_ids()
##                    candidate_links_tree = {}
##                    for candidate_link_id in candidate_links:
##                        candidate_links_tree[candidate_link_id] = []
##                    break
##                # Move back another spot.
##                else:
##                    self.note('[find_routepoint] Last routepoint has no link ID: backing up.')
##                    last_routepoint_id -= 1
##                    last_link_ids = self.routepoints[last_routepoint_id].get('link_ids', [])[:1]
##            # Contra argument would be if you've got no previous routepoint. (?)
##            if len(last_link_ids) > 0:
##                self.note('[find_routepoint] Last routepoint has link %i.' % last_link_ids[0])
##                candidate_links_tree = link_network.get_near_link_tree(
##                    link_id = last_link_ids[0],
##                    max_node_distance = cursor_size
##                    )
##                candidate_links = candidate_links_tree.keys()
##
##        self.note('[find_routepoint] Candidate links %s.' % str(candidate_links))
##        self.note('[find_routepoint] Candidate links tree %s.' % str(candidate_links_tree))
##
##        # Update the potential link vertices.
##        candidate_vertices = link_network.get_link_vertex_tuples(candidate_links)
##
##        # Find nearest candidate vertex.
##
##        # Trackpoint geometry.
##        xtp = self.trackpoints[trackpoint_id].get('x')
##        ytp = self.trackpoints[trackpoint_id].get('y')
##        ztp = self.trackpoints[trackpoint_id].get('z')
##
##        # Init near evaluators.
##        near_link_id = None
##        near_vertex_position = None
##        near_xyz = None
##        near_distance = None
##
##        for vertex in candidate_vertices:
##            xcv, ycv, zcv = vertex[2]
##            vertex_distance = ((xcv - xtp) ** 2 + (ycv - ytp) ** 2) ** 0.5
##            # If first checked, automatically gets assigned.
##            # If this vertex is nearer than already measured, assign.
##            if near_distance is None or vertex_distance < near_distance:
##                near_link_id = vertex[0]
##                near_vertex_position = vertex[1]
##                near_xyz = (xcv, ycv, zcv)
##                near_distance = vertex_distance
##
##        # List of link IDs that the routepoint represents should omit the last
##        # ID in the walk-back, since that's the last RP's link.
##        rp_link_ids = [near_link_id] + candidate_links_tree[near_link_id][:-1]
##
##        self.note('[find_routepoint] Routepoint links %s' % str(rp_link_ids))
##        self.note('[find_routepoint] Link vertex %i' % near_vertex_position)
##
##        return rp_link_ids, near_vertex_position, near_xyz
##
##
##    # Works (2012-05-21).
##    def get_od_node_ids(self):
##        """Returns the origin/destination node pair IDs as a tuple.
##
##        Requires that the entry & exit nodes be already determined. Currently,
##        entry & exit nodes are only added during the fix_routelink_backtracks
##        function, so be sure that is run first.
##
##        """
##
##        # No way of determining O&D nodes if there's no links!
##        if len(self.routelinks) == 0:
##            onode_id, dnode_id = None, None
##        else:
##            last_position = sorted(self.routelinks)[-1]
##            try:
##                onode = self.routelinks[1]['entry_node']
##                onode_id = onode[0]
##            # Means no entry node found on first link.
##            except KeyError:
##                onode_id = None
##            try:
##                dnode = self.routelinks[last_position]['exit_node']
##                dnode_id = dnode[0]
##            # Means no entry node found on final link.
##            except KeyError:
##                dnode_id = None
##
##        return onode_id, dnode_id
##
##
##    # Works (2012-05-21).
##    def get_routelink_at(self, position):
##        """Returns the link ID at the given position in the route."""
##
##        if position in self.routelinks:
##            result = self.routelinks[position].get('link_id')
##        else:
##            result = None
##
##        return result
##
##
##    # Works (2012-05-25).
##    def get_routepoint_at(self, position):
##        """Returns a routepoint's info."""
##
##        routepoint_info = {}
##
##        if self.routepoints.get(position):
##            for key in self.routepoints[position]:
##                routepoint_info[key] = self.routepoints[position].get(key)
##
##        return routepoint_info
##
##
##    # Works (2012-05-25).
##    def get_routepoint_ids(self):
##        """Returns all routepoint IDs."""
##
##        routepoint_id_list = []
##
##        for position in self.routepoints:
##            routepoint_id_list.append(position)
##
##        return routepoint_id_list
##
##
##    # Works (2012-05-21).
##    def get_routepoints(self):
##        """Returns all routepoints along the route."""
##
##        routepoint_list = []
##
##        for position in self.routepoints:
##            routepoint_list.append(self.routepoints.get(position))
##
##        return routepoint_list
##
##
##    # Works (2012-05-21).
##    def get_ordered_route_links(self):
##        """Returns an ordered list of link IDs."""
##
##        ordered_route_links = []
##
##        for link_position in sorted(set(position for position in self.routelinks)):
##            ordered_route_links.append(
##                self.routelinks[link_position].get('link_id')
##                )
##
##        return ordered_route_links
##
##
##    # Routes: Data management.
##
##    # Works (2012-07-06): Added multi-link routepoints.
##    def add_routepoint_to_route(self, trackpoint_id, link_ids,
##                                vertex_position, xyz):
##        """Writes link, vertex position & location to given trackpoint.
##
##        This effectively writes to a new trackpoint collection, one which moves
##        the points to run along the network.
##
##        """
##        self.routepoints[trackpoint_id] = {
##            'link_ids': link_ids, 'vertex_position': vertex_position,
##            'x': xyz[0], 'y': xyz[1], 'z': xyz[2]
##            }
##
##        return trackpoint_id
##
##
##    # Works (2012-06-29): Added multi-link routepoints.
##    def set_entire_route(self, link_network, avoid_problems=True,
##                         cursor_size=1):
##        """Sets all the routepoints & routelinks using a cursor-like method.
##
##        Cursor size limits the network nodal distance the cursor can check from
##        the previous routepoint's location.
##
##        The avoid_problems switch allows the cursor to skip trackpoints that
##        may throw the cursor off course. Problem flags are assigned via other
##        functions.
##
##        """
##        self.note('[set_entire_route] Start.')
##
##        # Initialize the routepoints.
##        self.routepoints = {}
##        # Add position 1, but with no values in it (needed for later get()).
##        self.routepoints[1]= {}
##
##        for trackpoint_id in sorted(self.trackpoints):
##
##            if avoid_problems:
##                # Has problems: Skip.
##                if self.trackpoints[trackpoint_id].get('problems', False):
##                    do_this_one = False
##                # Has no problems: Skip.
##                else:
##                    do_this_one = True
##            # Not avoiding, who cares?
##            else:
##                do_this_one = True
##
##            if do_this_one:
##                link_ids, vertex_position, xyz = self.find_routepoint(
##                    trackpoint_id, link_network, cursor_size
##                    )
##                self.add_routepoint_to_route(trackpoint_id, link_ids,
##                                             vertex_position, xyz)
##
##        # Align the routelinks with the newly-set routepoints.
##        self.set_routelinks()
##
##        # Fix backtrack problems in routelinks.
##        ##self.fix_routelink_backtracks(link_network)
##
##        self.note('[set_entire_route] End.')
##
##        return [self.routelinks[link]['link_id'] for link in self.routelinks]
##
##
##    # Works (2012-07-06): Added multi-link routepoints.
##    def set_routelinks(self):
##        """Creates the routelink dictionary.
##
##        This can only be effectively run after all routepoints are added to
##        the route.
##
##        """
##
##        route_links = []
##        link_position = 1
##        last_link_id = None
##
##        for routepoint_id in self.routepoints:
##            routepoint = self.routepoints.get(routepoint_id)
##            # Only add trackpoint info that were calculated into the route.
##            if routepoint:
##                self.note('[set_routelinks] Routepoint %i.' % routepoint_id)
##                # Add link_ids in in travel order (reverse the walk-back).
##                for this_link_id in reversed(routepoint['link_ids']):
##                    # Only add a link if preceding trackpoint wasn't there also.
##                    if this_link_id != last_link_id:
##                        route_links.append((link_position, this_link_id))
##                        link_position += 1
##                        self.note('[set_routelinks] Added link %i.' %
##                                  this_link_id)
##                    else:
##                        self.note('[set_routelinks] Skipping link %i.' %
##                                  this_link_id)
##                    last_link_id = this_link_id
##
##        # Order the tuples for link positon.
##        sorted_route_links = sorted(route_links)
##
##        for link_position, link_id in sorted_route_links:
##            self.routelinks[link_position] = {'link_id': link_id}
##
##        return
##
##
##    # Works (2012-07-09): Added multi-link routepoints; cleaned OD node code.
##    def fix_routelink_backtracks(self, link_network):
##        """Finds & corrects route links that cause a 'backtrack'.
##
##        If a routepoint is assigned to an incorrect side-link, it will appear
##        as if it has 'backtracked' from going down the side-link. In truth,
##        the side-link in question was never accessed, it's just a mix-up in the
##        assignment. This function will identify these mix-ups, and correct them.
##
##        This function will also add an entry_node & exit_node to routelinks.
##
##        Method of identification: When the route changes links, it cannot enter
##        from the same node it entered the previous link with.
##
##        Prescribed fix: Remove side-link from the route links, and re-evaluate
##        previous link.
##
##        """
##        self.note('[fix_routelink_backtracks] Start.')
##        messages = []
##        state = True
##
##        # Skip backtracks for no-link: Give O&D nodes NoneType tuples.
##        if len(self.routelinks) == 0:
##            self.note('[fix_routelink_backtracks] No-link trip.')
##            messages.append('No-link trip.')
##            self.routelinks[1]['entry_node'] = (None, None)
##            self.routelinks[1]['exit_node'] = (None, None)
##
##        # Skip backtracks for single-link: Assign O&D with T&F, respectively.
##        elif len(self.routelinks) == 1:
##            self.note('[fix_routelink_backtracks] Single-link trip.')
##            messages.append('Single-link trip')
##            fnode1, tnode1 = link_network.get_link_nodes(self.get_routelink_at(1))
##            self.routelinks[1]['entry_node'] = fnode1
##            self.routelinks[1]['exit_node'] = tnode1
##
##        # For all other link-counts, correct backtracks & assign true O&Ds.
##        else:
##
##            # Start with the first position.
##            this_position = 1
##            last_position = this_position - 1
##            next_position = this_position + 1
##
##            # Loop will run until the current link position doesn't exist.
##            while this_position in self.routelinks:
##
##                self.note('[fix_routelink_backtracks] Link position %i.' % (this_position))
##
##                # Get values for current & next link.
##                this_link_id = self.get_routelink_at(this_position)
##                if last_position in self.routelinks:
##                    last_link_id = self.get_routelink_at(last_position)
##                else:
##                    last_link_id = None
##                if next_position in self.routelinks:
##                    next_link_id = self.get_routelink_at(next_position)
##                else:
##                    next_link_id = None
##
##                fnode, tnode = link_network.get_link_nodes(this_link_id)
##
##                # Assign entry & exit nodes.
##                # Entry node on from-end.
##                if last_link_id in link_network.get_link_ids_at_node(fnode[0], fnode[1]):
##                    self.routelinks[this_position]['entry_node'] = fnode
##                # Entry node on to-end.
##                elif last_link_id in link_network.get_link_ids_at_node(tnode[0], tnode[1]):
##                    self.routelinks[this_position]['entry_node'] = tnode
##                # last_link_id is None (first position).
##                else:
##                    self.routelinks[this_position]['entry_node'] = (None, None)
##                # Exit node on from-end.
##                if next_link_id in link_network.get_link_ids_at_node(fnode[0], fnode[1]):
##                    self.routelinks[this_position]['exit_node'] = fnode
##                # Exit node on to-end.
##                elif next_link_id in link_network.get_link_ids_at_node(tnode[0], tnode[1]):
##                    self.routelinks[this_position]['exit_node'] = tnode
##                # next_link_id is None (last link position).
##                else:
##                    self.routelinks[this_position]['exit_node'] = (None, None)
##
##                this_entry_node = self.routelinks[this_position]['entry_node']
##                this_exit_node = self.routelinks[this_position]['exit_node']
##                if last_position in self.routelinks:
##                    last_exit_node = self.routelinks[last_position].get('exit_node')
##                else:
##                    last_exit_node = (None, None)
##
##                # Check for backtrack.
##
##                # Exit node is the same: This means we've backtracked.
##                # Or, last link_id is the same as this one.
##                # Ignore if last routelink.
##                if (this_exit_node == last_exit_node or
##                    this_link_id == last_link_id):
##                    # Remove the backtracked link.
##                    del self.routelinks[this_position]
##                    self.note('[fix_routelink_backtracks] Removed backtrack position %i link %i.' % (this_position, this_link_id))
##                    # Reset the routelink order to fill the gap.
##                    self.reset_routelink_positions()
##                    # Start process over from previous location.
##                    this_position -= 1
##                # Exit node different: This means a good transition, no backtrack.
##                else:
##                    # Positive result means we move to the next position.
##                    this_position += 1
##
##                last_position = this_position - 1
##                next_position = this_position + 1
##
##            # Set origin node.
##            try:
##                first_position = 1
##                first_id = self.get_routelink_at(first_position)
##                for node in link_network.get_link_nodes(first_id):
##                    if node != self.routelinks[first_position]['exit_node']:
##                        self.routelinks[first_position]['entry_node'] = node
##            except KeyError:
##                self.note('[fix_routelink_backtracks] KeyError while setting o-node.')
##                messages.append('KeyError while setting o-node')
##                state = False
##                try:
##                    self.routelinks[first_position]['entry_node'] = (None, None)
##                except KeyError:
##                    self.note('[fix_routelink_backtracks] KeyError: Backtracks removed all links.')
##                    messages.append('KeyError: Backtracks removed all links')
##
##            # Set destination node.
##            try:
##                final_position = [pos for pos in self.routelinks][-1]
##                final_id = self.get_routelink_at(final_position)
##                for node in link_network.get_link_nodes(final_id):
##                    if node != self.routelinks[final_position]['entry_node']:
##                        self.routelinks[final_position]['exit_node'] = node
##            except IndexError:
##                self.note('[fix_routelink_backtracks] IndexError: Backtracks removed all links.')
##                messages.append('IndexError: Backtracks removed all links')
##                state = False
##            except KeyError:
##                self.note('[fix_routelink_backtracks] KeyError while setting d-node.')
##                messages.append('KeyError while setting d-node')
##                state = False
##
##        self.note('[fix_routelink_backtracks] End.')
##
##        return state, messages
##
##
##    # Works (2012-05-21).
##    def reset_routelink_positions(self):
##        """Removes gaps in the link positions (1, 1 + 1, 1 + 2, ...)."""
##
##        old_route = self.routelinks
##        new_route = {}
##        new_link_position = 1
##
##        for old_link_position in sorted(position for position in old_route):
##
##            new_route[new_link_position] = old_route[old_link_position]
##            self.routelinks = new_route
##            new_link_position += 1
##
##        return
##
##
### Works (2012-05-21).
##def load_gps_trips_from_mysql(data_dictionary, trip_ids=None,
##                              logfile_path=None, trip_limit=None):
##    """Loads a given track from an MySQL database.
##
##    Requires the MySQLdb library to be installed, and a conversion
##    dictionary for the source data. Returns a list of the GPSTrip class
##    instances for the given trip IDs (optional: omitting gets all in DB) and
##    omitting all trips noted in the logfile (optional).
##
##    """
##    import csv
##    try:
##        import MySQLdb
##    except ImportError:
##        raise Exception
##
##    table_name = data_dictionary.get('source_table')
##    trip_id_field = data_dictionary.get('source_id_field')
##    timestamp_field = data_dictionary.get('source_timestamp_field')
##    easting_field = data_dictionary.get('source_easting_field')
##    northing_field = data_dictionary.get('source_northing_field')
##    elevation_field = data_dictionary.get('source_elevation_field')
##    xy_accuracy_field = data_dictionary.get('source_xy_accuracy_field')
##    z_accuracy_field = data_dictionary.get('source_z_accuracy_field')
##
##    # Connect to DB.
##    conn = MySQLdb.Connect(host = data_dictionary.get('source_host'),
##                           port = data_dictionary.get('source_port'),
##                           user = data_dictionary.get('source_user'),
##                           passwd = data_dictionary.get('source_password'),
##                           db = data_dictionary.get('source_db'),
##                           compress = 1)
##
##    # Gather all wanted trip IDs from the DB.
##    trips_from_database = []
##    cursor = conn.cursor(cursorclass = MySQLdb.cursors.DictCursor)
##    cursor.execute("select distinct %s from %s" % (trip_id_field, table_name))
##    gps_trip_ids = cursor.fetchall()
##    for row in gps_trip_ids:
##        # Including all trips.
##        if trip_ids is None:
##            trips_from_database.append(int(row[trip_id_field]))
##        # Including just the given trips.
##        else:
##            if int(row[trip_id_field]) in trip_ids:
##                trips_from_database.append(int(row[trip_id_field]))
##
##    trips_to_analyze = []
##
##    # Check the logfile (if included).
##    if logfile_path:
##        # Gather all logged trip IDs.
##        trips_logged = []
##        log_rows = csv.DictReader(open(logfile_path, 'rU'))
##        for row in log_rows:
##            trip_logged_as = row.get('analysis_complete', str(None))
##            if trip_logged_as in (str(True), str(False)):
##                trips_logged.append(int(row.get(trip_id_field)))
##        # Omit logged trips.
##        for trip_id in trips_from_database:
##            if trip_id in trips_logged:
##                pass
##            else:
##                trips_to_analyze.append(trip_id)
##    else:
##        trips_to_analyze = trips_from_database
##
##    # If a trip limit is set, truncate trips_to_analyze.
##    if trip_limit:
##        trips_to_analyze = trips_to_analyze[:trip_limit]
##
##    # Gather all rows for the chosen trips.
##    all_trips_rows = []
##    for trip_id in  trips_to_analyze:
##        cursor.execute(
##            "select * from %s where trip_id = %s order by %s, %s" %
##            (table_name, str(trip_id), trip_id_field, timestamp_field)
##            )
##        trip_rows = cursor.fetchall()
##        all_trips_rows.append(trip_rows)
##
##    # Close MySQL connection.
##    cursor.close()
##    conn.close()
##
##    # Create class instances and write trackpoints into it.
##    all_trips = []
##    for trip_rows in all_trips_rows:
##        if trip_rows:
##            data_dictionary['trip_id'] = trip_rows[0].get(trip_id_field)
##            this_trip = GPSTrip(data_dictionary)
##            track_position = 1
##            # Add each trackpoint to the trip class instance.
##            for row in trip_rows:
##                # Convert GPS lat/lon to state plane feet.
##                x, y = transform_gps_to_orsps(
##                    longitude = row.get(easting_field),
##                    latitude = row.get(northing_field)
##                    )
##                this_trip.add_trackpoint(
##                    position = track_position,
##                    timestamp = row.get(timestamp_field),
##                    xyz = (x, y, row.get(elevation_field)),
##                    accuracy = (row.get(xy_accuracy_field),
##                                row.get(z_accuracy_field))
##                    )
##                track_position += 1
##
##            all_trips.append(this_trip)
##
##    return all_trips