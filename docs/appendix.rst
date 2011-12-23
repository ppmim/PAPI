Appendix
========

.. _superalign:

``superalign``
**************

Algorithm
---------

``superalign`` operates in a very tree-like fashion:

1. First, the program groups the exposures in terms of those that substantially overlap ('GroupExpoures').


2. Then, group by group (pointing by pointing), it attempts to iteratively align all the exposures 
   (``AlignEntirePointing``).  As it successfully finds alignments that work between exposures within a 
   group, it generates a stacked star list -- where CRs are eliminated and contains an ever increasing  
   number of pointlike objects (individual exposures often contain more contamination and miss real objects).


3. After generating a best-effort stack of all the stars from each group, it begins to build up a 
   mosaic by starting at the central group and iteratively accreting nearby groups (to produce an ever 
   larger list of stars).  As each group is accreted, the relative position/rotation of all the groups 
   added to that point (within the ever growing mosaic) is  reoptimized to minimize the overall error.


``superalign`` uses a variant on the similar triangle method to
determine the alignment between two individual catalogs (see
``TryAlignment`` routine in simplealign.c).  However, instead of
attempting to find similar triangles in both images, ``superalign``
attempts to find similar bi-directional vectors (object pairs with the
same separation and the same relative orientation).  After finding
similar vectors, ``superalign`` uses these vectors to determine a
transformation from one catalog to the other.  Each transformation is
given a score based upon how well it maps objects onto each other.
The transformation with the highest score is then adopted.  One final
set of iterations (see ``FindSolution`` routine in simplealign.c) is
then attempted (starting with the transformation which had the highest
score) to try to further improve the score.

Operation
---------


After compiling with a "make superalign" command, ``superalign``
is currently invoked from the command line as follows:

::

    superalign INPUT_SHIFTS OUTPUT_STARLIST OUTPUT_SHIFTS 

The first file is the INPUT_SHIFTS file and should have the following format:

::

    [# of exposures]
    [exposure #1]  [x-shift guess #1] [y-shift guess #1] [rot-angle guess #1]
    [exposure #2]  [x-shift guess #2] [y-shift guess #2] [rot-angle guess #2]
    [exposure #3]  [x-shift guess #3] [y-shift guess #3] [rot-angle guess #3]
    ....

The first exposure in the INPUT_SHIFTS file gives the reference image,
and will define the final coordinate system.  It will have an x,y
shift equal to 0,0 and a rotation angle equal to 0.

The second file is the OUTPUT_STARLIST file and contains a compilation
of all the objects in the final mosaic on a common grid.  Following
the construction of this object list, it should be trivial to perform
alignments onto the common frame.  The format for the file is as follows:

::

    1 [x position #1] [y position #1] [magnitude #1]
    2 [x position #2] [y position #2] [magnitude #2]
    3 [x position #3] [y position #3] [magnitude #3]
    .....

The third file is the OUTPUT_SHIFTS file and a list of all the determined
shifts:

::

    [exposure #1]  [x-shift #1] [y-shift #1] [rot-angle #1]
    [exposure #2]  [x-shift #2] [y-shift #2] [rot-angle #2]
    [exposure #3]  [x-shift #3] [y-shift #3] [rot-angle #3]
    ....

For each input catalog, ``superalign`` will output a cleaned catalog.
The cleaned catalog (with the same file name as the input catalog
except with ".name" appended) will contain those objects ``superalign``
believes are real (from the stacking analysis it does at each
position).

Examples
--------

To demonstrate the usage of this package, an example has been provided
in the "examples" directory.  This directory contains a list of
(distortion free) catalogs.  Each catalog is contained in a ".matchin"
file and corresponds to a single ACS exposure which will contribute to
some final set of drizzled images.

For this example, the usage is:

::

    superalign supermatch.in supermatch.stars supermatch.out 

"supermatch.stars" will contain the final object list for the mosaic on
the common grid.

"supermatch.out" will contain the derived shifts and rotations for all
of the individual exposures.

Free Parameters
---------------

To improve performance, one can change the following parameters which
are part of the ``SuperAlignParmRec`` structure (change with caution!):

::

    float TypDist;      /* Typical RMS error for objects which are well aligned */
    float MaxDist;      /* MaxDist is equal to the distance two objects can be apart to even register as matches. Even if objects register as matches -- this does not imply they will have a huge effect on the final score. */
    float TolPair;      /* Objects with Neighbors Closer than TolPair pixels are thrown out (due to possible deblending problems) */
    float MaxShiftErr;  /* Maximum Pixel Error in X/Y Shifts */
    float MaxRotErr;    /* Maximum Error in Rotation Angle (Degrees) */
    float GroupSize;    /* Maximum Distance allowed between two pointings to be considered part of the same grouping */

X and Y (Distance) Units
------------------------

Arbitrary units can be used for the x/y positions and distances
throughout, but they must be consistent.  Among the shifts, x/y
positions which should be consistent are:

    * X/Y positions of objects in the catalogs (for the individual exposures)
    * Best guess x/y shifts given for each exposure in INPUT_SHIFTS file
    * The Five Following Parameters from SuperAlignParmProp (at the end of superalign.c):
        * TypDist
        * MaxDist
        * TolPair
        * MaxShiftErr
        * GroupSize
