"""Control script for conducting map-matching analysis on the CycleLane data.

This script is specifically configured for use with the GPS data collected from
the CycleLane mobile app's recorded trips. The trips are stored in a central
MySQL database. This script loads the GPS trackpoints directly from the
database into a custom Python class format. In this format, map-matching
analysis is performed, accompanied by some clean-up and data quality checking.

After analysis, this script writes-off a number of output formats, including
GIS data types for visualization and CSV-format text files for modeling
applications & statistical analysis.

Note: The message handler & functions imported from RouteModelingLib labelled
with '_arcgis'require the local computer to have the library arcpy installed.
This library is installed as part of Esri ArcGIS 10.x or greater.

Note: Functions imported from RouteModelingLib labelled with '_mysql'require
the local computer to have the library MySQLdb installed. The best place to get
the MySQLdb installer is at <codegood.com/archives/4>.

Note: Functions imported from RouteModelingLib labelled with '_from_mysql'
require the local computer to have the library pyproj installed. The best place
to get the pyproj installer is at <code.google.com/p/pyproj>.

"""

import os
import sys
import arcpy

sys.path.append(r'\\clsrv111\GIS\common\Python\modules')
import batchutilities as bu

sys.path.append(r'\\clsrv111\transpor\Models\BikeModels\LCOG\Scripts')
from RouteModelingLib_Breakin import (create_csv_od_from_arcgis,
                                      create_csv_routes_from_arcgis,
                                      create_output_layers_for_arcgis,
                                      get_all_trip_ids,
                                      load_gps_trips_from_mysql,
                                      load_link_network_from_arcgis,
                                      update_od_nodes_in_arcgis,
                                      update_routepoints_in_arcgis,
                                      update_routes_in_arcgis,
                                      update_trackpoints_in_arcgis)

# Model network version.
NET_VERSION = '2012'
# Combined version (use this if you're testing and need more specific runs).
VERSION = NET_VERSION
# Set to True if you want to use the GPSTrips internal debugger logging.
# Beware using the debugger logging for large amounts of trips, since it will
# exponentially increase the amount of memory used by a trip and disk writing.
LOG_TRIPS_DEBUGGING = False
# Analysis area envelope ((easting range), (northing range)).
# Using TMA boundary's outer limits right now.
AREA_ENVELOPE = ((4203973.0, 4304304.0), (4235928.0, 4255359.0))
# Distance of largest gap allowed, or smallest trip allowed.
DISTANCE_THRESHOLD = 528.0  # or 1/10 of a mile.
# Number of nodes the cursor can traverse for route assigning.
CURSOR_SIZE = 2
# Allows us to proportion trip-crunching for each script-run, and how much each
# 'chunk' loaded has. I've never had a single chunk greater than 512. I've also
# noticed 32-bit memory ceiling issues (1GB) with 256.
SCRIPT_TRIP_LIMIT = 512
CHUNK_TRIP_LIMIT = 64

# Workspaces.
DATA_PATH = r'\\clsrv111\transpor\Models\BikeModels\LCOG\Data'
VERSION_PATH = os.path.join(DATA_PATH, 'CycleLaneOutput_%s' % VERSION)
LAYER_PATH = os.path.join(VERSION_PATH, 'TripGISLayers_%s' % NET_VERSION)
MODEL_DATA_PATH = os.path.join(DATA_PATH, 'BikeModelData_%s.gdb' % NET_VERSION)
TRIP_PATH = os.path.join(VERSION_PATH, 'CycleLaneTrips_%s.gdb' % NET_VERSION)

# Files & classes.
LINK_NETWORK_PATH = os.path.join(MODEL_DATA_PATH, 'BikeFacility')
LOGFILE_PATH = os.path.join(VERSION_PATH,
                            'CycleLaneLogfile_%s.csv' % NET_VERSION)
OD_PATH = os.path.join(TRIP_PATH, 'OriginsDestinations')
OD_CSV_PATH = os.path.join(VERSION_PATH,
                           'CycleLaneOriginsDestinations_%s.csv' % NET_VERSION)
ROUTES_PATH = os.path.join(TRIP_PATH, 'Routes')
ROUTES_CSV_PATH = os.path.join(VERSION_PATH,
                               'CycleLaneRoutes_%s.csv' % NET_VERSION)
ROUTEPOINTS_PATH = os.path.join(TRIP_PATH, 'Routepoints')
TRACKPOINTS_PATH = os.path.join(TRIP_PATH, 'Trackpoints')

# Info dictionaries.
# Data loader info.
NETWORK_LOAD_INFO = {
    'source_path': LINK_NETWORK_PATH,
    'source_data_type': 'ArcGIS',
    'source_id_field': 'link_id',
    'source_from_node_field': 'fnode', 'source_to_node_field': 'tnode',
    'source_from_level_field': 'z_F', 'source_to_level_field': 'z_T',
    'shape_field': 'Shape',
    'source_subselect_sql': "bik_rest <> 1",
    }
TRIPS_LOAD_INFO = {
    'source_host': 'clsrv300.ris.lane.or.us',
    'source_port': 3306,
    'source_user': 'jb',
    'source_password': 'bikedb2011',
    'source_db': 'cycletracks',
    'source_table': 'coord',
    'source_data_type': 'MySQL',
    'source_id_field': 'trip_id',
    'source_timestamp_field': 'recorded',
    'source_easting_field': 'longitude',
    'source_northing_field': 'latitude',
    'source_elevation_field': 'altitude',
    'source_xy_accuracy_field': 'hAccuracy',
    'source_z_accuracy_field': 'vAccuracy'
    }
# ArcGIS output updaters.
LAYERS_CREATE_INFO = {
    'folder_path': LAYER_PATH,
    'name_format': 'Trip <ID>',  # No .lyr; added in script.
    'template_path': os.path.join(LAYER_PATH, 'Trip Template.lyr'),
    'sublayers_sql': {
        'OriginsDestinations': 'trip_id = <ID>',
        'Routes': 'trip_id = <ID>',
        'Routepoints': 'trip_id = <ID>',
        'Trackpoints': 'trip_id = <ID>'
        }
    }
OD_UPDATE_INFO = {
    'output_path': OD_PATH,
    'trip_id_field': ('trip_id',),
    'node_id_field': ('node_id',),
    'node_type_field': ('node_type',),
    'shape_field': 'Shape'
    }
ROUTEPOINTS_UPDATE_INFO = {
    'output_path': ROUTEPOINTS_PATH,
    'route_id_field': ('trip_id',),
    'routepoint_id_field': ('routepoint',),
    'link_ids_field': ('link_ids',),
    'link_vertex_field': ('link_vertex',),
    'shape_field': 'Shape'
    }
ROUTES_UPDATE_INFO = {
    'output_path': ROUTES_PATH,
    'route_id_field': ('trip_id',),
    'link_id_field': ('link_id',),
    'link_position_field': ('link_position',),
    'shape_field': 'Shape'
    }
TRACKPOINTS_UPDATE_INFO = {
    'output_path': TRACKPOINTS_PATH,
    'track_id_field': ('trip_id',),
    'trackpoint_id_field': ('trackpoint',),
    'link_ids_field': ('link_ids',),
    'link_vertex_field': ('link_vertex',),
    'shape_field': 'Shape'
    }
# CSV output updaters.
OD_CSV_UPDATE_INFO = {
    'output_path': OD_CSV_PATH,
    'trip_id_field': ('trip_id',),
    'node_id_field': ('node_id',),
    'node_type_field': ('node_type',),
    }
ROUTES_CSV_UPDATE_INFO = {
    'output_path': ROUTES_CSV_PATH,
    'route_id_field': ('trip_id',),
    'link_id_field': ('link_id',),
    'link_position_field': ('link_position',)
    }


# Works (2012-06-22).
def update_message(status_msgs, data_location):
    """Simple execution of data update status-messages."""

    if status_msgs[0]:
        mh.send('Updated %s at %s.' %
                (data_location, bu.getcurrenttime()))
    else:
        mh.send('%s update failed at %s.' %
                (data_location, bu.getcurrenttime()))
    if len(status_msgs[1]) > 0:
        mh.send('Messages: %s.' % '; '.join(status_msgs[1]))

    return


# Works (2012-06-22).
def trip_message(status, trip_id, trip_notes=[]):
    """Simple execution of trip analysis state messages."""

    if status:
        mh.send('Trip %i matched at %s.' % (trip_id, bu.getcurrenttime()))
    else:
        mh.send('Trip %i match failed at %s.' % (trip_id, bu.getcurrenttime()))
    if len(trip_notes) > 0:
        note_row = '%s.' % '; '.join(trip_notes)
        mh.send(note_row)
    else:
        note_row = ''
    mh.sendlog('%s,%s,%s' % (trip_id, status, note_row))

    return


# Works (2012-06-22).
def sort_logfile(logfile_path):
    """Performs a sort on the logfile. Sort fields are hardwired."""

    import csv
    from operator import itemgetter

    with open(logfile_path, 'rU') as logfile:
        logreader = csv.DictReader(logfile)
        logrows = sorted(logreader,
                         key = itemgetter('trip_id', 'analysis_complete', 'notes'))
    # Add field-name row.
    logrows.insert(0, {'trip_id': 'trip_id',
                       'analysis_complete':'analysis_complete',
                       'notes':'notes'})
    with open(logfile_path, 'wb') as logfile:
        logwriter = csv.DictWriter(
            logfile, fieldnames = ('trip_id', 'analysis_complete', 'notes')
            )
        logwriter.writerows(logrows)

    return


# Works (2012-06-29).
def main():
    """Controls the procedural order of the script when run by itself."""

    # Loading data.

    loaded_network = load_link_network_from_arcgis(NETWORK_LOAD_INFO)
    mh.send('Link network loaded version %s at %s.' % (NET_VERSION,
                                                       bu.getcurrenttime()))

    if loaded_network:

        trips_so_far = 0
        while trips_so_far < SCRIPT_TRIP_LIMIT:

            loaded_trips = load_gps_trips_from_mysql(
                data_dictionary = TRIPS_LOAD_INFO,
                # Specific trips only for testing.
                trip_ids = None,
                logfile_path = LOGFILE_PATH,
                trip_limit = CHUNK_TRIP_LIMIT
                )

            # Add these trips to the script count.
            trips_so_far += len(loaded_trips)


            # Ensures logfile lines will start on new row.
            mh.sendlog('')
            mh.send('%i GPS trips loaded from database at %s.' %
                    (len(loaded_trips), bu.getcurrenttime()))

            # Have some trips: Analyze each.
            if loaded_trips:

                for trip in loaded_trips:

                    trip_messages = []
                    is_good = True

                    while is_good:

                        # PRE-VET: Must be within the envelope of the area of
                        # study. Doesn't currently work: Need to work on
                        # envelope logic.
                        ##trackpoints_out = trip.find_trackpoints_outside_area(
                        ##    easting_range = AREA_ENVELOPE[0],
                        ##    northing_range = AREA_ENVELOPE[1]
                        ##    )
                        ##number_out = len(trackpoints_out)
                        ##if number_out > 0:
                        ##    trip_messages.append('PREFAIL-OUT')
                        ##    trip_messages.append('%i out-of-area: %s' %
                        ##                         (number_out,
                        ##                          str(trackpoints_out)))
                        ##    is_good = False
                        ##    break

                        # PRE-VET: Must be longer than the distance threshold.
                        tdistance, ttime = trip.find_absolute_trip_length()
                        if tdistance < DISTANCE_THRESHOLD:
                            trip_messages.append('PREFAIL-SHORT')
                            trip_messages.append(
                                'Total length of %0g feet in %i seconds' %
                                (tdistance, ttime)
                                )
                            is_good = False
                            break

                        # PRE-VET: Must not have gaps larger than the distance
                        # threshold.
                        gdistance, gtime = trip.find_largest_trackpoint_gap()
                        if gdistance > DISTANCE_THRESHOLD:
                            trip_messages.append('PREFAIL-GAP')
                            trip_messages.append(
                                'Gap of %2g feet over %i seconds' %
                                (gdistance, gtime)
                                )
                            is_good = False
                            break

                        # Set route & fix backtrack issues.
                        trip.set_entire_route(link_network = loaded_network,
                                              avoid_problems = True,
                                              cursor_size = CURSOR_SIZE)
                        is_good, messages = trip.fix_routelink_backtracks(
                            link_network = loaded_network
                            )
                        trip_messages.extend(messages)

                        break

                    trip_message(is_good, trip.trip_id, trip_messages)

                    if LOG_TRIPS_DEBUGGING:
                        with open(os.path.join(VERSION_PATH,'DebugTrip%i.txt' %
                                               trip.trip_id), 'wb') as dfile:
                            for trip_note in trip.all_notes:
                                dfile.write('%s\n' % trip_note)

                # Write output data (updates).

                # Trackpoints (ArcGIS).
                status_msgs = update_trackpoints_in_arcgis(
                    gps_trips = loaded_trips,
                    output_dictionary = TRACKPOINTS_UPDATE_INFO
                    )
                update_message(status_msgs, TRACKPOINTS_PATH)

                # Routes (ArcGIS).
                status_msgs = update_routes_in_arcgis(
                    gps_trips = loaded_trips,
                    link_network = loaded_network,
                    output_dictionary = ROUTES_UPDATE_INFO
                    )
                update_message(status_msgs, ROUTES_PATH)

                # Routepoints (ArcGIS).
                status_msgs = update_routepoints_in_arcgis(
                    gps_trips = loaded_trips,
                    output_dictionary = ROUTEPOINTS_UPDATE_INFO
                    )
                update_message(status_msgs, ROUTEPOINTS_PATH)

                # O&D nodes (ArcGIS).
                status_msgs = update_od_nodes_in_arcgis(
                    gps_trips = loaded_trips,
                    link_network = loaded_network,
                    output_dictionary = OD_UPDATE_INFO
                    )
                update_message(status_msgs, OD_PATH)

                # Layer files (ArcGIS).
                # Use the idlist below to only (re)create these trips' layers.
                idlist = [trip.trip_id for trip in loaded_trips]
                # Use the idlist below to (re)create all trips' layers.
                ##idlist = get_all_trip_ids(TRACKPOINTS_PATH, 'trip_id')
                status_msgs = create_output_layers_for_arcgis(
                    trip_ids = idlist,
                    output_dictionary = LAYERS_CREATE_INFO
                    )
                update_message(status_msgs, LAYER_PATH)

                # Need to close & open log so next chunk can read updated rows.
                mh.closelog()
                mh.filehandle = open(LOGFILE_PATH, 'a')

            # Have no more trips loaded: Break out of loop.
            else:
                break

        # Logfile cleanup.

        mh.closelog()
        mh.send('All output located at: %s' % VERSION_PATH)
        sort_logfile(LOGFILE_PATH)
        mh.send('Logfile: %s' % os.path.basename(LOGFILE_PATH))

        # Write output data (whole rewrites).

        # O&D nodes (CSV).
        status_msgs = create_csv_od_from_arcgis(
            source_path = OD_PATH,
            data_dictionary = OD_CSV_UPDATE_INFO
            )
        update_message(status_msgs, OD_CSV_PATH)

        # Routes (CSV).
        status_msgs = create_csv_routes_from_arcgis(
            source_path = ROUTES_PATH,
            data_dictionary = ROUTES_CSV_UPDATE_INFO
            )
        update_message(status_msgs, ROUTES_CSV_PATH)


# Works (2012-06-22).
if __name__ == '__main__':

    try:
        # Construct the message handler.
        mh = bu.messagehandler(log = True,
                               filename = LOGFILE_PATH,
                               append = True)
        # Console message setting.
        mh.echo = True
        # Mail message settings.
        mh.mail = False
        ##mh.sender = 'jblair@lcog.org'
        ##mh.recipient = 'jblair@lcog.org'
        ##mh.subject = '%s Log for %s' % (bu.program, bu.gettoday())

        # Logfile only for recording the individual trip messages.
        # Turn off auto-write using this, use sendlog(message) for logfile.
        mh.logfile = False

        mh.send('%s started at %s.' % (bu.program, bu.getcurrenttime()))

        main()

        status = 0

    except:
        # Handle errors by setting the exit status (so if a batch file wraps
        # this script can use it). A utility function formats the messages from
        # the geoprocessor.
        status = 1
        bu.handleexception(arcpy, mh)

    finally:
        # Cleanup actions include sending a mail message to the administrator
        # to indicate program success or failure, closing the logfile, and
        # returning the exit status.

        mh.send('%s %s at %s.' % (bu.program, bu.StatusMessage[status],
                                  bu.getcurrenttime()))

    mh.send('Program status: %i (0=good, 1=bad)' % status)
    mh.closelog()
    raw_input('Press Enter key to exit')
    ##sys.exit(status)
