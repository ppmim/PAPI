.. _os:

Operating System
****************

Installation
============
Once you have configured the computer system, the next step is to install the operating system `openSUSE 13.1 <https://en.opensuse.org/Portal:13.1>`_.
To do this, you only have to follow the instructions found `here <https://en.opensuse.org/SDB:DVD_installation>`_.

Basically, the steps are:

Boot
----
Insert the openSUSE13.1 x64 DVD into the DVD drive and restart the computer. After the power-on self test (POST), 
the boot screen should appear. Use the keyboard up and down buttons to select Installation and press enter to confirm.

.. note::

    The power-on self test (POST) can last several minutes, mainly because the kernel must recognize the PowerEdge RAID Controller (PERC) H710 card.
    

.. image:: _static/131-installer-1.png
   :align: center
   :scale: 65 %

  

The kernel will load and some debugging messages may appear. These can usually be ignored unless there is an issue 
starting the installation. In that event, these messages may make debugging easier. 

.. image:: _static/131-installer-2.png
   :align: center
   :scale: 65 %


Welcome
-------

This page provides a choice of languages and keyboard layouts. The installer will choose a keyboard layout automatically after selecting a language. 
For PANIC we will select English (US) as Language and Spanish(ES) as Keyboard Layout.

The license agreement is below the language options. Language translations are available below. Click Next to accept the license and continue. 

.. image:: _static/131-installer-3.png
   :align: center
   :scale: 65 %
   
Installation Mode
-----------------

This guide only covers new installations. Leave the default options and click Next to continue.

Clock and Time Zone
-------------------

Set the computer's time zone to **UTC**, which can be selected from the drop down box below the map.

Leave Hardware Clock Set to UTC enabled.

If the Date and Time shown below the menu is incorrect, click Change to enter the correct date and UTC time manually.

If the computer has an active internet connection, select Synchronize with NTP Server to obtain a 
date and time automatically from a time server. The default server is usually the best selection. Later we will set the CAHA NTP server (derfel.caha.es). 
If the Date and Time shown below the menu is incorrect, click Change to enter the correct date and time manually.

Desktop Selection
-----------------

Select `KDE Desktop` and click Next.


Partitioning
------------

By default, the installer will provide a suggested partition setup, however we will select ``Edit Partition Setup``. 
It will open the `Expert Partitioner` where drive layouts and installation can be customized. The partitions to 
create are::

    - /dev/sda1             swap (   2 GB)
    - /dev/sda2   /         ext3 ( 456 GB)
    - /dev/sdb1   /data1    ext3 ( 1,8 TB)
    - /dev/sdc1   /data2    ext4 ( 5,5 TB)
    
    
For swap partition, we use the default suggested size of only 2 GB (in principle we do not need a larger swap parittion 
because we have a large amount of RAM and we are not going to use the suspend mode to allow the contents of the 
computer's RAM to be stored in disk).


Create New Users
----------------
Create next users::
    
    - panic     /home/panic   (bash shell) --> for tests
    - obs22     /home/obs22   (bash shell) --> for observations
    

and in both cases, disable ``Automatic Login`` to be prompted for a password when the computer starts. 

Installation Settings
---------------------
An installation summary will appear before installation begins. Confirm the options selected and click Install to begin the process. Clicking on each of the underlined titles will bring up the respective settings page.


131-installer-13.png 

Software Selection
------------------

Selecting the Software title will bring up the Software Selection screen, where groups of packages (known as patterns) can be selected. For a desktop installation, the default options are usually acceptable. Software packages can always be added or removed after the installation.

Select Details to manage all available packages individually, rather than as patterns. 
