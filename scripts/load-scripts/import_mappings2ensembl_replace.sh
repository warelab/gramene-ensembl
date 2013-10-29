#!/bin/bash

REG='/usr/local/gramene/conf/ensembl.registry'
SCRIPT='/usr/local/gramene/scripts/markers/import_mapping_markersdb2ensembl.pl'
#TRACK_SP='Oryza => Rice'


SP=$1
MAPSET=$2
TRACK_SP=$3

if test ! "$SP"
then
    echo "Need a species in the registry ${REG}"
    exit;
fi

if test ! "$MAPSET"
then
    echo "Need a markers db mapset name"
    exit;
fi

if test ! $TRACK_SP
then
    echo "Need a track_sp such as 'Oryza => Rice' "
    exit;
fi
echo "DATE IS $DATE; SPECIES is "$SP" (in ${REG}); Markersdb type is $TYPE; markers db mapset is "$MAPSET"; track_species is "${TRACK_SP}"";


perl $SCRIPT -registry_file $REG -species "$SP" -mapset_name "$MAPSET" -marker_type "EST Cluster"  -cs toplevel -replace

perl $SCRIPT -registry_file $REG -species "$SP" -mapset_name "$MAPSET" -marker_type RFLP -track_species "${TRACK_SP}"  -cs toplevel -replace

perl $SCRIPT -registry_file $REG -species "$SP" -mapset_name "$MAPSET" -marker_type OVERGO  -track_species "${TRACK_SP}" -cs toplevel -replace

perl $SCRIPT -registry_file $REG -species "$SP" -mapset_name "$MAPSET" -marker_type Microarray_Probe  -track_species "${TRACK_SP}" -cs toplevel -replace

