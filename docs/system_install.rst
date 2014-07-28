
.. _system_installation:

Installation & Configuration  
============================

This document describes the procedure to install the PANIC computer from scratch, 
and shows how to install the PANIC software.

This guide is intended for technical audiences and CAHA personnel working in 
the field of IT.

Main steps:

- Raid configuration
- OS installation
- System Packages installation
- GEIRS installation
- OT installation
- PAPI installation


System overview
***************

PANIC computer (PC) consist in a Dell PowerEdge R720 two-socket 2U 
rack servers with the next key components:

* Processor:
  Intel(R) Xeon(R) CPU E5-2643 0 @ 3.30GHz (QPI 8,0GT/s, 10 M caché , 130W, 4C)

* Memory:
  8 x 8GB RDIMM, 1600 MHz, 

* Disks:
  - 2 x 500GB SATA 7200 3,5" HD hot-plug
  - 6 x 2 TB SAS nearline 6 Gb/s 7,2K rpm 3,5" HD hot-plug

* BIOS: 
  Dell 1.6.0

* Video: 
  Matrox Graphics, Inc. G200eR2

* Network: 
  Broadcom 5720 Quad-port Gigabit Ethernet

* RAID Controller
  PERC H710p integrated , 1GB NV Cache

* Power Supply:
  Dual, Hot-plug, Redundant Power Supply (1+1), 1100W

* Optical device:
  16X DVD+/-RW SATA

* iDRAC7 Enterprise card

* 7 PCIe slots: 
    - One x16 full-length, full-height 
    - Three x8 full-length, full-height 
    - Three x8 half-length, half-height



3Yr Basic Warranty - Next business day
5Yr ProSupport Next Business day On-site service

* Service Tag: 5KNBWX1
* Express Service Code: 12132422245


In addition, the system will be installed with the home made (MPIA) OPT-PCIe 
board, the hardware interface for the data acquisition from the Read-Out 
electronics. 


Note: Due to PANIC can be installed both 2.2m telescope and 3.5m telescope,
we have two indentical systems (Dell R720), one for each telescope. This way,
each system can be used as a spare system for the other, with only a slight 
change on the software setup (GEIRS, OT and PAPI).


RAID configuration
******************
The system has 8 hard disks, 2x500GB for OS in a RAID-1 configuration, 2x2TB 
for GEIRS data and 4x2TB for the data processing pipeline. The Raid configuration
for the syste is described below:


Virtual Disk 1 - RAID 1 (OS system +  PANIC software) – 500 MB
Disk 0  - 500 GB SATA
Disk 1  - 500 GB SATA
Stripe Element Size - 64KB (default)
Write Policy - Write-Back (default)
Read Policy - Adaptative-Read-Ahead (default)

Virtual Disk 2 - RAID 1 (Data acquisition directory) – 2 TB
Disk 2  - 2 TB Nearline SAS
Disk 3  - 2 TB Nearline SAS
Stripe Element Size - 1 MB (non default)
Write Policy - Write-Through  (non default)
Read Policy - Adaptative-Read-Ahead (default)

Virtual Disk 3 - RAID 5  (Data pipeline products) - 6TB
Disk 4  -  2 TB Nearline SAS
Disk 5  -  2 TB Nearline SAS
Disk 6  -  2 TB Nearline SAS
Disk 7  -  2 TB Nearline SAS

Stripe Element Size - 1 MB (default)
Write Policy - Write-Back (default)
Read Policy - Adaptative-Read-Ahead (default)


Virtutal disk creation procedure
--------------------------------

1 Start or restart the Dell PowerEdge server.
2 Press Ctrl+R to access the PERC 710p Integrated BIOS Configuration Utility.
The configuration utility opens and shows the VD Mgmt tab.
3 To create a new virtual disk, press F2 and select Create New VD.
4 Select an appropriate RAID configuration from the RAID Level drop-down menu.
In our case, we use RAID1 for VD 1 and 2, and RAID5 for virtual disk 3.
5 Under Physical Disks, select the hard disks that you want to include in the virtual disk.
6 Select OK and press Enter.
The virtual disk is created.
7 Expand the Disk Group option, and under Virtual Disks, select the newly created virtual disk and press F2 to open the Operations menu.
8 From the Operations menu, select Initialization > Start Init. to initialize the new virtual disk. It will take several hours.
The new RAID volume is ready for use.


RAID-1   465.25 GB   HDD Adaptive Read Ahead Write Back  64K No  
RAID-1   1862.50 GB  HDD Adaptive Read Ahead Write Through   1M  No  
RAID-10  3725.00 GB  HDD Adaptive Read Ahead Write Back  1M  No  






What is the iDRAC?
******************
The iDRAC (integrated Dell Remote Access Controller) is an integrated tool which
can perform a multitude of configuration and maintenance tasks. By logging on 
to the iDRAC, you can fix many problems yourself immediately, instead of having 
to rely on our support staff.

How do I connect to the iDRAC?

There are two IP addresses; one of them for the server itself and the other one 
for the iDRAC.
[...]



OPT-PCIe board
**************

Installing OS
*************

1 If the server is currently running, power down the server.

2 Insert the installation DVD that you burned from the openSuSE 13.1 ISO image 
into the server’s DVD drive.

3 Power up the server.
The server boots from the DVD and displays the OS type selection screen.
It can takes a long time, so be patient and wait for the display of the initial 
presentation menu.

4 Select Installation

5 Select your language, keyboard layout and accept the licence terms:
    
    Language: English (US)
    Keyboard: Spanish

6 Installation Mode: we choose a new installation

7 The installer analyzes your hardware and builds the software repository cache:

8 Clock and Time Zone

Set the timezone here. It's recommended to set the hardware clock to UTC. 

9 Desktop Selection
Select KDE

10 Partitioning

By default openSUSE will propose to create three new partitions / (root) for 
system files, /home/ for personal files of users and swap which is used as a 
supplement for RAM, 

11 Create New User
Now it's time to create PANIC user. Note that by default the root user 
(administrator) password will be the same as the password for the PANIC user.

If you want the added security of a separate root password, consider unchecking 
that checkbox. You may also want to consider disabling autologin to prevent 
people from easily accessing your system and data.

12 Installation Settings 

- Double check that everything is as desired - this is the point of no return!

- Scroll down to the Firewall and SSH section and enable SSH.

- Check under the Software section that at least are selected the next software 
packages:

TODO

13 Actual Installation

Now the actual installation is performed.

14 Automatic Configuration

After installation is performed, the system will restart and perform 
autoconfiguration. And finally your brand new openSUSE system will start. 



References
**********
- PowerEdgeand R720xd and R720 Technical Guide: 
http://i.dell.com/sites/doccontent/shared-content/data-sheets/en/Documents/dell-poweredge-r720-r720xd-technical-guide.pdf
- Dell PowerEdge RAID Controller (PERC) H710P User’s Guide
- Dell PowerEdge 12th Generation Server BIOS Configuration
- Updating BIOS on DELL 12G PowerEdge Servers
- openSUSE Start-Up Guide: http://activedoc.opensuse.org/book/opensuse-start-up
- openSUSE Reference Guide: http://activedoc.opensuse.org/book/opensuse-reference
- openSUSE System Analysis and Tuning Guide: http://activedoc.opensuse.org/book/opensuse-system-analysis-and-tuning-guide
- openSUSE Security Guide: http://activedoc.opensuse.org/book/opensuse-security-guide
- R. J. Mathar, PANIC - Generic Infrared Software - Graphical User Interface, 
PANIC-SW-DCS-01.pdf URL: http://www.mpia.de/~mathar/public/PANIC-SW-DCS-01.pdf 

.. _PANIC: http://www.iaa.es/PANIC
.. _CAHA: http://www.caha.es
.. _Omega2000: http://www.caha.es/CAHA/Instruments/O2000/index.html
.. _HAWK-I: http://www.eso.org/sci/facilities/paranal/instruments/hawki/
.. _sphinx: http://sphinx.pocoo.org
.. _pdf: http://www.iaa.es/~jmiguel/PANIC/PAPI/PAPI.pdf
  