

def map_num_to_letter(df, col='Row', inplace=True):
    """
    Maps integer row numbers (1–26) to uppercase alphabet letters.

    Args:
        df (pd.DataFrame): Input DataFrame.
        col (str): Column to remap. Defaults to 'Row'.
        inplace (bool): Modify df in place (True) or return a copy (False).

    Returns:
        pd.DataFrame 
    """

    import string

    row_map = {
        i: letter for i, 
               letter in enumerate(string.ascii_uppercase, start=1)
    }

    mapped = df[col].map(row_map) # Map the column values using the row_map

    if inplace:
        df[col] = mapped
        return None
    else:
        result = df.copy()
        result[col] = mapped
        return result


class Categorizer:

    """Make a threshold for low and high EV-uptakers. 
    Next, Label the samples based on their function. The threshold is 
    determined using the untreated (Control) group."""
    
    def __init__(self, normalized_column:str, df_untreated=None, 
                 df_treated=None, func="mean", sd=2, 
                 lower_outlier="Inhibitors", normal="Normal", 
                 higher_outlier="Inducers"):
        
        """ 
        Args: 
            df_untreated: The control group
            df_treated: The sample group 
            normalized_column: The column that contains the normalized values 
                of the variable that you want to used for defining the 
                threshold.
            func: The aggregation function for calculating the threshold. 
                It may be "mean" or "median".
            sd: The number of standard deviation away from the mean/median. 
                You may use 1, 2, 3 for 1SD, 2SD, and 3SD, respectively. 
            lower_outlier: The value(name) you want to give to the group of 
                values below the lower threshold. The default value is 
                "Inhibitors".
            normal: The value(name) you want to give to the group of values 
                bewtween the two lower and upper threshold level. 
                The default value is "Normal".
            higher_outlier: The value(name) you want to give to the group of 
                values above the upper threshold. 
                The default value is "Inducers".

        Returns:
            pd.DataFrame: The treated dataframe including a new 
                category column based on the SD threshold.
        """
        import pandas as pd 
    
        self.df_untreated = df_untreated
        self.df_treated = df_treated
        self.normalized_column = normalized_column
        self.func = func
        self.sd = sd
        self.lower_outlier = lower_outlier
        self.normal = normal
        self.higher_outlier = higher_outlier
        
        
    def threshold_generator(self):
        
        # This is either mean or median of untreated (control), 
        # depending on the desired parameter
        self.total_parameter_control = self.df_untreated[
            self.normalized_column
        ].agg(self.func) 
        
        # Make 1sd or sd, depending on the input of sd in the function
        self.std_control = (
            self.df_untreated[self.normalized_column].std()
        ) * self.sd 
        
        self.plus_sd = (
            self.total_parameter_control + self.std_control
        ) # upper threshold
        self.minus_sd = (
            self.total_parameter_control - self.std_control
        ) # lower threshold
        
        print(
            f"The values for control {self.func}, {self.func}-{self.sd}SD "
            f"and {self.func}+{self.sd}SD, respectively:  "
            f"{self.total_parameter_control}, {self.minus_sd}, "
            f"{self.plus_sd}"
        )
        
        # the control-related values: mean, lower_thresh, higher_thresh
        return self.total_parameter_control, self.minus_sd, self.plus_sd 
    
    def category_generator(self):
        self.df_treated = self.df_treated.copy()
        cat_col = 'Category' + "_by_" + str(self.sd) + "SD"
        
        self.df_treated.loc[
            self.df_treated[self.normalized_column] < self.minus_sd, 
            cat_col
        ] = self.lower_outlier
        
        self.df_treated.loc[
            (self.df_treated[self.normalized_column] >= self.minus_sd) & 
            (self.df_treated[self.normalized_column] <= self.plus_sd), 
            cat_col
        ] = self.normal
        
        self.df_treated.loc[
            self.df_treated[self.normalized_column] > self.plus_sd, 
            cat_col
        ] = self.higher_outlier

        inhibitors = len(
            self.df_treated.loc[self.df_treated[cat_col] == self.lower_outlier]
        )
        inducers = len(
            self.df_treated.loc[self.df_treated[cat_col] == self.higher_outlier]
        )
        total = len(self.df_treated)

        print(
            (f"***Based on {self.sd}SD*** \nFrom the total rows number of "
             f"{total:,} in treated group: \n{inhibitors:,} or"),
            (f"{(inhibitors/total)*100:.2f}% of the rows are inhibitors \n"
             f"and \n{inducers:,} or {(inducers/total)*100:.2f}% of the "
             f"rows are inducers"),
            f"\n\nTHE NEW DATAFRAME IS GENERATED."
        )
        return self.df_treated


# Caclulate the slope of lines for gating 
def line_formula(x1, y1, x2, y2):
    """Args:
        x1, y1: Coordinates of the first point.
        x2, y2: Coordinates of the second point.
    Returns:
        m: Slope of the line.
        b: y-intercept of the line.
    """
    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1
    print("The line equation is: y = {:.2f}x + {:.2f}".format(m, b))
    return m, b # return slope and y-intercept

def position(row, m_line, b_line, x_axis_thresh, y_axis_thresh,
            x_axis_col= 'Nuclei - Intensity Nucleus Alexa 488 Mean',
            y_axis_col= 'Nuclei - Intensity Nucleus Alexa 568 Mean'): 

    """
    Args:
        row: A row of a DataFrame containing the y-coordinate and the 
            y-coordinate on the line.
        y_coord_col: The name of the column containing the actual y-coordinate 
            of the point.
        y_on_line_col: The name of the column containing the y-coordinate 
            calculated from the line formula using the actual x-coordinate.
    Returns:
        The position of data point in relation to a line.
    """
    if row[x_axis_col] < x_axis_thresh and row[y_axis_col] < y_axis_thresh:
        return 'Below_threshold'
    else:
        # Calculate the y-coordinate on the line using the x-coordinate of 
        # the point.
        y_line = (m_line * row[x_axis_col]) + b_line 
        
        # Compare the y-coordinate of the point with the y-coordinate on the 
        # line
        if row[y_axis_col] > y_line:
            return 'above'
        elif row[y_axis_col] < y_line:
            return 'below'
        else:
            # Either simply 'on' or you can use np.isclose() or to check if 
            # two numbers are close enough to be considered equal
            # Or row[y_coord_col] == row[y_on_line_col]
            return 'on' 

def cytotox_group(row, position1_col, position2_col):
    """
    Args:
        row: each row of the dataframe containing the cytotoxicity value. 
            This will be passed to apply(), later.
        position1_col: The col with position of the point compared to the 
            first (left) line.
        position2_col: The col with position of the point compared to the 
            second (right) line.
    Returns:
        Determines which group of cytotixicity (apoptosis, late apoptosis, 
            necrosis) the datapoint belongs to.
    """
    if (row[position1_col] == 'Below_threshold' and 
            row[position2_col] == 'Below_threshold'):
        return 'Viable'
    
    else:
        if row[position1_col] == 'above' and row[position2_col] == 'above':
            return 'Necrosis'
        elif row[position1_col] == 'below' and row[position2_col] == 'above':
            return 'Late_Apoptosis'
        elif row[position1_col] == 'below' and row[position2_col] == 'below':
            return 'Apoptosis'
        else:
            return 'Unknown'



def cytotox_scatter_plot(df, casp_threshold, pi_threshold, 
                         x1_line1, y1_line1, x2_line1, y2_line1,
                         x1_line2, y1_line2, x2_line2, y2_line2, 
                         incubation = "30min",
                         combine_duplicates=True, 
                         width=1300, height=850, 
                         facet_col_spacing=0.02,
                         facet_row_spacing=0.09,
                         x_axis_tick_label_font_size=14,
                         y_axis_tick_label_font_size=14,
                         x_axis_label_x_position=0.5,
                         x_axis_label_y_position=-0.1,
                         y_axis_label_x_position=0.45,  
                         y_axis_label_y_position=-0.07, 
                         x_axis_label_font_size=20,
                         y_axis_label_font_size=20,
                         sep_line_x_axis_position=3.5,  
                         sep_line_start_y_position=0.0,
                         sep_line_end_y_position=30000,
                         title = False,
                         ): 
    """A function for generating scatterplot for the result of cytotoxicity 
    assay.

    Args:
        df: The dataframe as the input.
        casp_threshold: The threshold for caspase 3/7.
        pi_threshold: The threshold for PI.
        x1_line1, y1_line1, x2_line1, y2_line1: The coordinates of the 
            first line.
        x1_line2, y1_line2, x2_line2, y2_line2: The coordinates of the 
            second line.
        incubation: The incubation time of the cells. You can choose between 
            "30min" and "16h" to change the graph title accordingly.
        combine_duplicates: Whether to combine the duplicates in the same plate
            format or not. No special calculation will be done for the 
            duplicates before combining.
    Returns:
        fig (plotly.graph_objects.Figure): A scatterplot for the result of 
            cytotoxicity assay.
    """


    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    if combine_duplicates:
        facet_col = "Sample"
    else:
        facet_col = "Plate format"


    fig = px.scatter(
    data_frame=df, x="Nuclei - Intensity Nucleus Alexa 488 Mean",
    y="Nuclei - Intensity Nucleus Alexa 568 Mean",
    facet_col=facet_col, facet_col_wrap=6,
    labels={
        "Nuclei - Intensity Nucleus Alexa 488 Mean":
        "Caspase 3/7 average signal intensity",
        "Nuclei - Intensity Nucleus Alexa 568 Mean":
        "PI mean signal intens."
    },
    color="cytotox_group",
    color_discrete_map={
        "Viable": "gray",  # or "black"
        "Apoptosis": "gray",  # or "blue"
        "Necrosis": "gray",  # "red"
        "Late_Apoptosis": "gray",  # or "#FF8C00"
    },
    facet_row_spacing=facet_row_spacing,
    facet_col_spacing=facet_col_spacing,
    width=width, height=height
)

    # Control the opacity (hue) and size of the spots on the graph.
    fig.update_traces(marker=dict(opacity=0.6, size=5))
    fig.update_xaxes(
        showline=True, linewidth=2, linecolor='black', mirror=True,
        ticks="outside", showticklabels=True, range=[0, 70000],
        tickfont=dict(size=x_axis_tick_label_font_size, color='black')
    )
    # fig.update_xaxes(ticksuffix = "  ") 
    # increase the distance of the ticklabel from the axis
    fig.update_yaxes(
        showline=True, linewidth=2, linecolor='black', mirror=True,
        ticks="outside", range=[0, 70000],
        tickfont=dict(size=y_axis_tick_label_font_size, color='black')
    )

    # Dash lines for showing the group of viable cells
    fig.add_shape(  # Adding the horizontal line for the threshold of PI
        type="line",
        x0=0, y0=pi_threshold,  # y_start: lower y-value
        x1=casp_threshold, y1=pi_threshold,    # y_end: upper y-value
        line=dict(color="blue", width=2),
        # dash="dot" for changing from default solid line to dash
        row="all", col="all"  # Apply to all subplots
    )

    fig.add_shape(  # Adding horizontal line for threshold of Caspase 3/7
        type="line",
        x0=casp_threshold, y0=0,  # y_start: lower y-value
        x1=casp_threshold, y1=pi_threshold,    # y_end: upper y-value
        line=dict(color="blue",  width=2),
        row="all", col="all"  # Apply to all subplots
    )

    # Add gating lines
    fig.add_shape(
        type='line',
        x0=x1_line1, y0=y1_line1,
        x1=x2_line1, y1=y2_line1,
        line=dict(color='red', width=1),  # dash='dash'
        row='all', col='all'
    )

    fig.add_shape(
        type='line',
        x0=x1_line2, y0=y1_line2,
        x1=x2_line2, y1=y2_line2,
        line=dict(color='green', width=1),  # dash='dash'
        row='all', col='all'
    )

    # Add bar to separate the subplots
    # Add a black vertical bar between the third and fourth subplot
    if sep_line_x_axis_position is not None:
        fig.add_shape(
            type="line",
            # Adjust this value to position the bar
            x0=sep_line_x_axis_position,
            y0=sep_line_start_y_position,    # Start at the bottom of y-axis
            x1=sep_line_x_axis_position,  # Same x-coordinate for vertical
            # End at top of y-axis (adjust based on y-axis range)
            y1=sep_line_end_y_position,
            line=dict(color="black", width=3),  # Black line with width of 3
            xref="paper",  # Use paper coordinates for x-axis
            yref="y domain"  # Use y-axis domain for y-axis
        )

    # change the order of traces drawn. 
    # This means that the blue spots (with lower data points than red spots) 
    # will be shown on top of the red spots. Otherwise, most of the of 
    # blue spots will remain hidden behind the crowd of red spots. 
    # Note1: Here, only changing the alpha value does not help. 
    # Note2: The order of tracer drawing is based on the order of values 
    # in the data. Here., the data of 1 h is normally drawn before 3 h. 
    # Nevertheless, we changed the drawing order.  
    fig.data = fig.data[::-1]

    # After changing the order of traces, the order of markers also 
    # changes in the legend. However, we can reverse the order back to the
    fig.update_layout(legend={'traceorder': 'reversed'}, )
    # Change background color of graphs and the paper color.
    fig.update_layout(plot_bgcolor="white", paper_bgcolor="white")
    # fig.update_layout(yaxis = dict(tickfont = dict(size=20)))
    fig.update_annotations(font_size=15)

    fig.add_annotation(
        x=casp_threshold + 2000, y=60000, text='Necrotic', font_size=10,
        textangle=0, showarrow=False, row="all", col="all"
    )
    fig.add_annotation(
        x=45000, y=pi_threshold + 4000, text='Apoptotic', font_size=10,
        textangle=0, showarrow=False, row="all", col="all"
    )
    fig.add_annotation(
        x=50000, y=60000, text='Late apoptotic', font_size=10,
        textangle=0, showarrow=False, row="all", col="all"
    )


    plate_to_sample = dict(zip(df["Plate format"], df["Sample"]))
    fig.for_each_annotation(
        lambda a: a.update(
            text=plate_to_sample.get(
                a.text.split("=")[-1].strip(),  # Extract value after "="
                a.text.split("=")[-1].strip()   # Fallback to just the value
            )
        )
    )

    # Remove axis titles from individual subplots
    fig.for_each_xaxis(lambda x: x.update(title=''))
    fig.for_each_yaxis(lambda y: y.update(title=''))

    # Add shared axis labels as annotations
    fig.add_annotation(
        text="Caspase 3/7 average signal intensity",  # Shared x-axis label
        xref="paper", yref="paper",
        x=x_axis_label_x_position, y=x_axis_label_y_position,
        showarrow=False,
        font=dict(size=x_axis_label_font_size, color="black"),
    )
    fig.add_annotation(
        text="PI mean signal intensity",  # Shared y-axis label
        xref="paper", yref="paper",
        x=y_axis_label_x_position, y=y_axis_label_y_position,
        showarrow=False,
        textangle=-90,
        font=dict(size=y_axis_label_font_size, color="black"),
    )
    if title:
        if incubation == "30min":
            title_text = "30-minute incubation"
        elif incubation == "16h":
            title_text = "16-hour incubation"
        else:
            title_text = "The title is not determined correctly"
        # Add a single title on top of the figure in center of the page.
        fig.update_layout(
            title_text=title_text,
            title_x=0.5,
            title_y=0.99,  # Adjust vertical position of the title
            title_font=dict(size=20, color="red"),
            # Add extra space at the bottom or top of the whole figure
            margin=dict(t=100)
        )

    fig.update_layout(showlegend=False)

    # Return the figure object for further use if needed.
    return fig  
    

def cytotox_count_stacked_bar(df, cytotox_column="cytotox_group", 
                              sample_column="Plate format", 
                              bar_width=0.3, combine_duplicates=True,
                              x_axis_tick_label_font_size=17, 
                              y_axis_tick_label_font_size=17,
                              x_axis_title_font_size=20, 
                              y_axis_title_font_size=20, title=False,
                              ):
    """
    Create a stacked bar graph of cytotoxicity groups by sample using raw 
    counts, with samples ordered by Samples_order column.

    Args:
        df (pandas.DataFrame): DataFrame containing cytotoxicity data.
        cytotox_column (str): Column name for cytotoxicity groups. 
            Defaults to "cytotox_group".
        sample_column (str): Column name for sample names. 
            Defaults to "Plate format".

    Returns:
        fig (plotly.graph_objects.Figure): Stacked bar graph figure.
    """
    import pandas as pd
    import plotly.graph_objects as go
    # Define cytotoxicity groups and colors
    cytotox_groups = [
        ("Late_Apoptosis", "#e4a11f"),
        ("Apoptosis", "green"),
        ("Necrosis", "#b30000"),
        ("Viable", "#23236d")
    ]
    
    # Create sample order mapping
    # Get unique samples using drop_duplicates as we only need one row 
    # for each well and ignore the rest of the cell measurements in that well. 
    # Then, from the unique rows, we only keep the columns that we need.
    sample_order_df = df.drop_duplicates(
        subset=[sample_column])[[sample_column, 'Sample', 'Samples_order']]

    # Sort the samples based on the Samples_order column, 
    # which already has our desired orders using the numbers starting from 1. 
    sample_order_df = sample_order_df.sort_values('Samples_order')
     
    # Get ordered lists
    # Get the ordered list of plate formats based on Samples_order column.
    sorted_plate_formats = sample_order_df[sample_column].tolist()
    # Get ordered list of samples based on the Samples_order column.
    sorted_samples = sample_order_df['Sample'].drop_duplicates().tolist()
    

    if combine_duplicates:
        # 1. Create a mapping from plate format to Sample
        plate_to_sample = dict(zip(df[sample_column], df["Sample"]))

        # 2. Group and reindex as before
        grouped_data = (
            df.groupby([sample_column, cytotox_column])
            .size()
            .reindex(
                pd.MultiIndex.from_product(
                    [sorted_plate_formats, [g[0] for g in cytotox_groups]],
                    names=[sample_column, cytotox_column]
                )
            )
            .fillna(0)
            .reset_index(name='count')
        )
        
        # 3. Replace plate format with values from Sample
        grouped_data['Sample'] = grouped_data[sample_column].map(
            plate_to_sample
        )

        # 4. Drop the original plate format column and reorder columns
        grouped_data = grouped_data.drop(columns=[sample_column])
        grouped_data = grouped_data[['Sample', cytotox_column, 'count']]

        # 5. Calculate mean of counts for duplicates (if any)
        grouped_data = (
            grouped_data
            .groupby(['Sample', cytotox_column], as_index=False)['count']
            .mean().rename(columns={'count': 'average_count'})
        )

        sample_order = dict(zip(df["Sample"], df["Samples_order"]))
        grouped_data['Samples_order'] = grouped_data['Sample'].map(
            sample_order
        )
        grouped_data = grouped_data.sort_values('Samples_order')
        grouped_data = grouped_data.drop(columns=['Samples_order'])

    else:
        # Group data while preserving order
        # This function creates a list of all possible combination of values 
        # of columns, ensuring all combinations are present in groupby df.
        grouped_data = (
            df.groupby([sample_column, cytotox_column])
            .size()
            .reindex(pd.MultiIndex.from_product(
                # The second parameter is the ordered list of groups.
                [sorted_plate_formats, [g[0] for g in cytotox_groups]],
                names=[sample_column, cytotox_column]
            ))
            .fillna(0)
            .reset_index(name='count')
        )

    # Create figure
    fig = go.Figure()
    
    # Add traces in reverse order for proper stacking
    for group, color in cytotox_groups:
        group_data = grouped_data[grouped_data[cytotox_column] == group]
        if combine_duplicates:
            x = group_data['Sample']
            y = group_data['average_count']
            yaxis_title = 'Average Count'
        else:
            x = group_data[sample_column]
            y = group_data['count']
            yaxis_title = 'Count'

        fig.add_trace(go.Bar(
            x=x,
            y=y,
            name=group,
            marker_color=color,
            hovertemplate="<b>%{text}</b><br>Count: %{y}<extra></extra>",
            width=bar_width,
        ))

    if title:
        title = 'Cytotoxicity Groups by Sample'
        fig.update_layout(title=title)
    
    # Update layout
    fig.update_layout(
        barmode='stack',
        xaxis_title='Sample',
        yaxis_title=yaxis_title,
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis={
            'tickvals': (sorted_samples if combine_duplicates 
                         else sorted_plate_formats),
            'ticktext': sorted_samples,
            'tickangle': 90,
            'type': 'category'  # Ensure categorical treatment
        },
        legend={'title': 'Cell Status'}
    )

    # Adding spines to x and y axes
    fig.update_xaxes(
        showline=True,      # Show axis line
        linewidth=1,        # Line width
        linecolor='black',  # Line color
        mirror=False,       # Do not mirror axis line on top
        showgrid=False,     # Hide gridlines
        zeroline=False,     # Hide zero line
        # Font size for x-axis labels
        tickfont=dict(color='black', size=x_axis_tick_label_font_size),
        # Font size for x-axis title
        title_font=dict(color='black', size=x_axis_title_font_size),
    )

    fig.update_yaxes(
        showline=True,
        linewidth=1,
        linecolor='black',
        mirror=False, # Do not mirror axis line on right
        showgrid=False,
        zeroline=False,
        tickfont=dict(color='black', size=y_axis_tick_label_font_size),
        # Font size for y-axis title
        title_font=dict(color='black', size=y_axis_title_font_size),
        range=[0, 2000]
    )

    return fig


def cytotox_percent_stacked_bar(df, cytotox_column="cytotox_group", 
                                sample_column="Plate format"):
    """
    Create a stacked bar graph of cytotoxicity groups by sample using 
    percentages, with samples ordered by Samples_order column.
    """
    import pandas as pd
    import plotly.graph_objects as go   
    # Define cytotoxicity groups and colors
    cytotox_groups = [
        ("Late_Apoptosis", "#e4a11f"),
        ("Apoptosis", "green"),
        ("Necrosis", "#b30000"),
        ("Viable", "#23236d")
    ]
    
    # Create sample order mapping
    sample_order_df = df.drop_duplicates(
        subset=[sample_column]
    )[[sample_column, 'Sample', 'Samples_order']]
    sample_order_df = sample_order_df.sort_values('Samples_order')
    
    # Get ordered lists
    sorted_plate_formats = sample_order_df[sample_column].tolist()
    sorted_samples = sample_order_df['Sample'].tolist()
    
    # Group data while preserving order
    grouped_data = (
        df.groupby([sample_column, cytotox_column])
        .size()
        .reindex(pd.MultiIndex.from_product(
            [sorted_plate_formats, [g[0] for g in cytotox_groups]],
            names=[sample_column, cytotox_column]
        ))
        .fillna(0)
        .reset_index(name='count')
    )

    # Calculate percentages
    total_counts = grouped_data.groupby(sample_column)['count'].transform('sum')
    grouped_data['percentage'] = grouped_data['count'] / total_counts * 100

    # Create figure
    fig = go.Figure()
    
    # Add traces in reverse order for proper stacking
    for group, color in cytotox_groups:
        group_data = grouped_data[grouped_data[cytotox_column] == group]
        
        fig.add_trace(go.Bar(
            x=group_data[sample_column],
            y=group_data['percentage'],
            name=group,
            marker_color=color,
            hovertemplate="<b>%{text}</b><br>Percentage: %{y:.2f}%<extra></extra>",
        ))
    
    # Update layout
    fig.update_layout(
        barmode='stack',
        title='Cytotoxicity Groups by Sample (Percentage)',
        xaxis_title='Sample',
        yaxis_title='Percentage',
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis={
            'tickvals': sorted_plate_formats,
            'ticktext': sorted_samples,
            'tickangle': 90
        },
        legend={'title': 'Cell Status'}
    )
    
    return fig