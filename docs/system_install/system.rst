.. _system:

System
******

PANIC computer system is compound by the **Dell PowerEdge R720** 12th Generation. It is a 2-socket, 2U server that features the 
Intel Xeon E5-2600 processor family and supports up to 768GB of DDR3 memory. Dell offers the R720 in various backplane configurations 
with up to 16 2.5-inch internal hard drives or 8 3.5-inch drives. 

The PANIC computer specifications are:

- 2 x Intel Xeon E5-2643 (3.3 GHz, QPI 8.0 GT/s, 10M cache, 130W, 4C) Turbo
- 3.5" chassis with up 8 Hard Drives
- 8 x 8 GB RDIMM 1600 MHz
- PERC H710p RAID Controller 1 GB cache
- 2 x 500 GB SATA 7200 3.5" 
- 6 x 2 TB SAS nearline 6 Gb/s 7.2K rpm 3.5"
- 16X DVD+/-RW unidad SATA
- iDRAC7 Enterprise
- 3Yr Basic Warranty - Next Business Day
- Dell Broadcom BCM5720 Quad-port 1 GB Network card
- 6 PCIe x8 + 1 PCIe x16 slots
- 5Yr ProSupport and Next Business Day On-Site Service


There are two identical computers, one for each telescope (2.2m and 3.5m). At the same time, each one can works as a spare of the other. The 
**service tag** of the computers are::
    
    - panic22 => BKNBWX1
    - panic35 => 5KNBWX1


.. image:: _static/Dell_PowerEdge_R720_inside.jpg
   :align: center
   :scale: 80 %

    

Installation
============

For the complete installation description of System, please go to 
`Dell PowerEdge R720 and R720xd Owner's Manual <http://www.dell.com/support/manuals/us/en/19/topic/poweredge-r720/720720XDOM-v3/en-us>`_.

RAID configuration
==================
The disk array configuration for PANIC computer is compound of three disk volumes:

- Disk volume 1 - RAID 1 (456 GB): for the Operating System and PANIC software binaries
- Disk volume 2 - RAID 1 (1.8 TB): for raw data (GEIRS)
- Disk volume 3 - RAID 5 (5.5 TB): for processed data from pipeline (PAPI)

Next figure describe that configuration:

.. image:: _static/panic_disks.png
   :align: center
   :scale: 80 %


Next steps describe the creation of the RAID(s) using the PERC H710 integrated BIOS configuration utility:

1. Launch PERC H710 Integrated BIOS Configuration Utility

After the disks are inserted, reboot the server. When the server is starting up, press Ctrl-R to 
launch the PowerEdge Expandable RAID Controller BIOS. Press Ctrl-R when it is displaying 
the following message on the console:

.. image:: _static/dell720_1-press-ctrl-r.png
   :align: center
   :scale: 80 %

This will launch the H710 Integrated BIOS Configuration Utility. This utility will have 
the following three TABs on the top.

- VD Mgmt – Virtual Disk Management, which is selected by default
- PD Mgmt – Physical Disk Management
- Ctrl Mgmt – Controller Management

You’ll see the 8 new disks that we added in the “Unconfigured Physical Disks” 
section under “VD Mgmt” tab.

  
2. Launch the Operations Menus

From the Virtual Disk Management, use arrow key and select 
‘PERC H710 Integrated (Bus 0x03, Dev 0x00)’. Press F2 to show available operations. 
This will display a pop-up menu with following choices. Select “Create New VD”:

- Create New VD
- Clear Config
- Foreign Config
- Manage Preserved Cache
- Security Key Management


3. Create a Virtual Disk with RAID-1 (**for Operating System and PANIC software**):

This will display a new screen. In it, do the following:

Press Enter on the RAID option, which will display the available RAID choices. In this case, I choose RAID-1.

Use Tab key and scroll to the “Physical Disks” section, which will display all the new disks available.

Press space-bar which will select the highlighted disks. In this case, we select the 1st two disks. The “X” 
in front of the disk indicates that it is selected to be included in this particular virtual disk creation.

In the “01:00:00″ and “01:00:01″ that is displayed in the front of the disks, the last two digits 
represent the slot number. This indicates that we have selected the disks that were inserted in 
the slot number 0 and 1 for this virtual disk creation.

Based on the number of disks selected, and the RAID type selected, the “VD Size” will show you how 
much usable disk space you will be getting out of this virtual disk. Since we used RAID-1, this will 
give usable disk space of only one disk (the other disk will be used for mirroring the data). 

Keep all the following options to the default values, and press enter when ‘OK’ is selected::

    - Element Size: 64KB
    - Read Policy: Adaptive Read
    - Write Policy: Write Back
    - Force WB with no battery unchecked
    - Initialize checked
    - Configure HotSpare unchecked

You might get a pop-up with the following message. Just click ‘OK’ on it:

“It is recommended that all newly created logical drives be initialized unless 
you are attempting to recreate a previous configuration and recover data as 
initialization is a destructive process.”


4. Review the new RAID-1 Disk Group

This will create the “Disk Group 1″ (with RAID-1). This displays the following information:

- Virtual Disks: This indicates that there is only one virutal disk in this disk group, with size of roughly 500 GB.
- Physical Disks: This indicates all the physical disks that are part of this virtual disk group.
- Total Free Capacity: This indicates the free capacity available in the virutal disks
- Hot Spare: This is empty in this example.

Since we have used 2 disks out of the 8 new disks that we inserted, we only see 6 disks under 
the “Unconfigured Physical Disks” section.

5. Create a Virtual Disk with RAID-1 (**for raw data from GEIRS**)

Use the arrow key and select “PERC H710 Integrated (Bus 0x03, Dev 0x00)”, Press F2, select “Create New VD”.

Again, use RAID-1 and select the next two disks (2 x 2 TB). Since we used RAID-1, this will give usable disk space 
of only one disk (the other disk will be used for mirroring the data). Keep all the following options 
to the next values, and press enter when ‘OK’ is selected::

    - Element Size: 1 MB
    - Read Policy: Adaptive Read
    - Write Policy: Write Through
    - Force WB with no battery unchecked
    - Initialize checked
    - Configure HotSpare unchecked

6. Review the new RAID-1 Disk Group

This will create the “Disk Group 2″ (with RAID-1). This displays the following information:

- Virtual Disks: This indicates that there is only one virutal disk in this disk group, with size of roughly 1.8 TB.
- Physical Disks: This indicates all the physical disks that are part of this virtual disk group.
- Total Free Capacity: This indicates the free capacity available in the virutal disks
- Hot Spare: This is empty in this example.

Since we have used 2 disks out of the 6 unused disks, we only see 4 disks under 
the “Unconfigured Physical Disks” section.

7. Create a Virtual Disk with RAID-5 (**for data from Pipeline and Quick-Look**)

Use the arrow key and select “PERC H710 Integrated (Bus 0x03, Dev 0x00)”, Press F2, select “Create New VD”.

This time, we use RAID-5 and select the next four disks (4 x 2 TB). Since we used RAID-5, this will 
give usable disk space of only 3 disks (the other disk will be used for data parity). 
Keep all the following options to the next values, and press enter when ‘OK’ is selected::

    - Element Size: 1 MB
    - Read Policy: Adaptive Read
    - Write Policy: Write Back
    - Force WB with no battery unchecked
    - Initialize checked
    - Configure HotSpare unchecked

8. Review the new RAID-5 Disk Group

This will create the “Disk Group 3″ (with RAID-5). This displays the following information:

- Virtual Disks: This indicates that there is only one virutal disk in this disk group, with size of roughly 5.5 TB.
- Physical Disks: This indicates all the physical disks that are part of this virtual disk group.
- Total Free Capacity: This indicates the free capacity available in the virutal disks
- Hot Spare: This is empty in this example.

Since we have used 4 disks out of the 4 unused disks, we do not see any other disks under 
the “Unconfigured Physical Disks” section.

9. Physical Disk Management (PD Mgmt)

To manage the physical disks themselves, Press Ctrl-N, which will take you to the “PD Mgmt” tab 
and display all the disks that are available on the server.

To perform any operations on these disks, use arrow keys and select a disk, and press “F2″. 
For most part, you don’t have to do anything here. Just be aware of what operations can be 
done on the disks, just in case, if you need it for any future use.



Once we have created the disk volumes, we are ready to install the Operating System.



