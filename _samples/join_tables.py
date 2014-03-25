def join(table1, table2):

    """
    NOT_USED : does not work, fails in atpy library; a bug ??
    """
    
    import numpy as np

    # Check that primary key is set in Table 1
    if table1._primary_key is None:
        raise Exception("Primary key of table 1 has to be set")

    # Check that primary key is set in Table 2
    if table2._primary_key is None:
        raise Exception("Primary key of table 2 has to be set")

    # Find all unique keys between the two tables
    keys = np.unique(np.hstack([table1[table1._primary_key],
                                table2[table2._primary_key]]))

    # Sort by increasing key
    keys.sort()

    # Find list of columns in final table (and meta-data)

    columns = []
    dtype = []
    units = []
    descriptions = []
    formats = []
    nulls = []

    for column in table1.columns:
        columns.append(column)
        dtype.append((column, table1.data[column].dtype))
        units.append(table1.columns[column].unit)
        descriptions.append(table1.columns[column].description)
        formats.append(table1.columns[column].format)
        nulls.append(table1.columns[column].null)

    for column in table2.columns:
        if column in columns:
            table2.rename_column(column, column + ".2")
            # raise Exception("Column %s already exists in Table 1" % column)

    for column in table2.columns:
        if column != table2._primary_key:
            columns.append(column)
            dtype.append((column, table2.data[column].dtype))
            units.append(table2.columns[column].unit)
            descriptions.append(table2.columns[column].description)
            formats.append(table2.columns[column].format)
            nulls.append(table2.columns[column].null)

    # Need to take into account vector columns
    # Need to use column header object, rather than doing it this long way

    dtype = np.dtype(dtype)

    # Create the new table and set it up
    table = atpy.Table()
    table._setup_table(len(keys), dtype=dtype, units=units,
                                  descriptions=descriptions,
                                  formats=formats, nulls=nulls)

    # Place the primary keys in the table
    table[table1._primary_key][:] = keys[:]
    table.set_primary_key(table1._primary_key)

    # Create lookup table to match IDs in Tables 1 and 2 to new Table
    t1_keys = dict(zip(table1[table1._primary_key], range(len(table1))))
    t2_keys = dict(zip(table2[table2._primary_key], range(len(table2))))

    # Find row number of Tables 1/2 containing rows of new Table (in order)
    from_table1 = [t1_keys[key] if key in t1_keys else -1
                   for key in table[table._primary_key]]
    from_table2 = [t2_keys[key] if key in t2_keys else -1
                   for key in table[table._primary_key]]

    # Convert to arrays
    from_table1 = np.array(from_table1)
    from_table2 = np.array(from_table2)

    # Figure out which rows in the new table have a match in Tables 1/2
    in_table1 = from_table1 >= 0
    in_table2 = from_table2 >= 0

    # Only use these to match data
    from_table1 = from_table1[in_table1]
    from_table2 = from_table2[in_table2]

    # Copy over the data to the new table
    for column in table1.columns:
        if column != table1._primary_key:
            table.data[column][in_table1] = table1.data[column][from_table1]
    for column in table2.columns:
        if column != table2._primary_key:
            table.data[column][in_table2] = table2.data[column][from_table2]

    return table


