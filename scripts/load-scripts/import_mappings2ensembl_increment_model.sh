#!/bin/bash

REG='/usr/local/gramene/conf/ensembl.registry'
SCRIPT='/usr/local/gramene/scripts/markers/import_mapping_markersdb2ensembl.pl'
#TRACK_SP='Oryza => Rice'


DATE=$1
SP=$2
MAPSET=$3
TRACK_SP=$4

if test ! $DATE
then
    echo "Need a date"
    exit;
fi

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

if test ! "$TRACK_SP"
then
    echo "Need a track_sp such as 'Oryza => Rice' "
    exit;
fi
echo "DATE IS $DATE; SPECIES is $SP (in ${REG}); Markersdb type is $TYPE; markers db mapset is $MAPSET; track_species is ${TRACK_SP}";


perl $SCRIPT -registry_file $REG -species "$SP" -mapset_name "$MAPSET" -marker_type "EST Cluster"  -track_species "Oryza => Rice" -date $DATE -cs toplevel

perl $SCRIPT -registry_file $REG -species "$SP" -mapset_name "$MAPSET" -marker_type "EST Cluster"  -track_species "Arabidopsis => Arabidopsis" -date $DATE -cs toplevel

perl $SCRIPT -registry_file $REG -species "$SP" -mapset_name "$MAPSET" -marker_type "EST Cluster"  -track_species "Zea => Maize" -date $DATE -cs toplevel

perl $SCRIPT -registry_file $REG -species "$SP" -mapset_name "$MAPSET" -marker_type "EST Cluster"  -track_species "monocot => monocot" -date $DATE -cs toplevel

perl $SCRIPT -registry_file $REG -species "$SP" -mapset_name "$MAPSET" -marker_type "EST Cluster"  -track_species "dicot => dicot" -date $DATE -cs toplevel


perl $SCRIPT -registry_file $REG -species "$SP" -mapset_name "$MAPSET" -marker_type RFLP  -track_species  "${TRACK_SP}" -date $DATE -cs toplevel

perl $SCRIPT -registry_file $REG -species "$SP" -mapset_name "$MAPSET" -marker_type OVERGO  -track_species "${TRACK_SP}" -date $DATE -cs toplevel

perl $SCRIPT -registry_file $REG -species "$SP" -mapset_name "$MAPSET" -marker_type Microarray_Probe  -track_species "${TRACK_SP}" -date $DATE -cs toplevel

