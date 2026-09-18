header <- dashboardHeader(

  titleWidth = 360,

title = tags$div(

  style = "
    display:flex;
    align-items:center;
    gap:10px;
    height:50px;
  ",

  tags$img(
    src = "logo.png",
    style = "
      height:49px;
      width:auto;
      flex-shrink:0;
    "
  ),

  tags$div(

    style = "
      display:flex;
      flex-direction:column;
      justify-content:center;
      align-items:flex-start;
      line-height:1;
    ",

    tags$div(

      HTML(
        "<span style='
            color:white;
            font-weight:700;
            font-size:26px;
          '>TPP</span>

         <span style='
            color:#0EA5A5;
            font-weight:300;
            font-size:26px;
            margin-left:1px;
          '>Analyst</span>"
      )

    ),
    tags$div(

style = "
width:100%;
height:1px;
background:#0EA5A5;
margin-top:3px;
margin-bottom:4px;
opacity:0.9;
"
),

tags$div(

  "THERMAL PROTEOMICS PLATFORM",

  style = '
      color:rgba(255,255,255,0.75);
      font-size:8px;
      letter-spacing:0.06em;
      margin-top:-1px;
      width:100%;
      text-align:left;
      white-space:nowrap;
    '

)

  )

)

)



sidebar <- dashboardSidebar(

  sidebarMenu(

    id = "tabs",
tags$li(
id = "sidebarHandle",
style = "
text-align:right;
padding:12px 14px;
color:white;
cursor:pointer;
list-style:none;
border-bottom:1px solid rgba(255,255,255,0.08);
",
icon("angles-left")
),

    menuItem("Welcome", tabName = "welcome", icon = icon("house")),
    menuItem("Data Input", tabName = "data_input", icon = icon("database")),
    menuItem("Quality Control", tabName = "QC", icon = icon("search")),
    menuItem("Normalization", tabName = "normalization", icon = icon("sliders")),
    menuItem("Clustering", tabName = "clustering", icon = icon("diagram-project")),
    menuItem("Model Fitting", tabName = "model_fitting", icon = icon("chart-simple")),
    menuItem("Protein Plots", tabName = "protein_plots", icon = icon("chart-column")),
    menuItem("GO terms", tabName = "go_terms", icon = icon("dna")),
    menuItem(
      "Comparison Models",
      tabName = "condition_models",
      icon = icon("scale-balanced")
    )
  )


  
)


body <- dashboardBody(

tags$head(

tags$script(HTML("

$(document).on('click', '#sidebarHandle', function () {

    $('body').toggleClass('sidebar-collapse');

    if ($('body').hasClass('sidebar-collapse')) {

        $('#sidebarHandle i')
            .removeClass('fa-angles-left')
            .removeClass('fa-angles-left fa-angles-right')
            .addClass('fa-angle-right');  

    } else {

        $('#sidebarHandle i')
            .removeClass('fa-angles-right')
            .removeClass('fa-angle-right')
            .addClass('fa-angles-left');
    }

});

")),


tags$link(
  rel = "stylesheet",
  href = "https://fonts.googleapis.com/css2?family=Geist:wght@300;400;500;600;700&display=swap"
),

  tags$style(HTML("

/* =========================================================
   TPP Analyst Theme Tokens
   Change colors ONLY here
   ========================================================= */

:root {

    /* -------- Branding -------- */

    --color-primary: #1E3A5F;
    --color-primary-hover: #17304E;

    --color-accent: #0EA5A5;
    --color-accent-hover: #0D8F8F;

--color-accent-soft: rgba(14,165,165,0.08);

   --color-secondary: #1E3A5F;

--color-accent-soft: rgba(14,165,165,0.10);


    --sidebar-icon: #9CA3AF;



    /* -------- Layout -------- */

    --sidebar-bg: #0F172A;
   --header-bg: #1E3A5F;

    --page-bg: #F7F7F8;

    /* -------- Cards -------- */

    --card-bg: #FFFFFF;
    --card-header-bg: #FFFFFF;

    --card-border: #ECECEC;
    --color-primary-ring: rgba(181,54,122,0.18);

    /* -------- Text -------- */

    --text-primary: #1E3A5F;
    --text-secondary: #6B7280;

    --text-light: #D1D5DB;

    /* -------- Status -------- */

    --success: #10B981;
    --warning: #F17C4A;
    --danger: #EF4444;

    /* -------- Tables -------- */

    --table-header: #1E3A5F;
--table-hover: #F5F6F8;

    --table-border: #F3E8FF;

    /* -------- Tabs -------- */

    --tab-bg: #F1F1F2;
    --tab-hover: #E8E8EA;

    /* -------- Tooltips -------- */

    --tooltip-bg: #1E3A5F;
    /* -------- UI States -------- */

--color-primary-ring: rgba(30,58,95,0.18);

--color-primary-soft: rgba(30,58,95,0.15);


/* -------- Typography -------- */

--font-family: 'Geist', sans-serif;

--tab-text: #4B5563;

--tooltip-highlight: #FBB44A;

/* -------- Borders -------- */

--input-border: #D1D5DB;


}

.sidebar-toggle {
    display: none !important;
}

.landing-hero {
  text-align: center;
  padding: 30px 40px 10px 40px;
}

.landing-hero h1 {
  font-size: 42px;
  font-weight: 700;
  color: #0f172a;
}

.landing-hero p {
  font-size: 18px;
  color: #64748b;
}
.landing-subtitle{

  margin-left:1px;     /* move right */

  margin-top:-1px;     /* move up/down */

  font-size:13px;      /* size */

  letter-spacing:0.08em; /* spacing between letters */


}
.landing-divider{

  width:100%;

  height:1px;

  background:#0EA5A5;

  margin-top:4px;

  margin-bottom:5px;

  opacity:0.9;

}
.workflow-card {
  background: white;
  border-radius: 16px;
  padding: 25px;
  text-align: center;
  box-shadow: 0 2px 8px rgba(228, 208, 208, 0.08);
}

.workflow {
  font-size: 22px;
  font-weight: 600;
}

.landing-card {
  background: white;
  border-radius: 16px;
  padding: 15px;
  box-shadow: 0 2px 8px rgba(0,0,0,.08);
  min-height: 340px;
}

.workflow-node{
  text-align:center;
  padding:20px;
}

.workflow-node i{
  font-size:20px;
  color:#0EA5A5;
  margin-bottom:8px;
}
.workflow-row{
  display:flex;
  align-items:center;
  justify-content:center;
  flex-wrap:wrap;
  gap:8px;
  margin-bottom:25px;
}

.workflow-node{
  display:flex;
  flex-direction:column;
  align-items:center;
  justify-content:center;

  width:75px;

  font-size:11px;
  text-align:center;
  color:#334155;
}

.workflow-node i{
  font-size:24px;
  margin-bottom:6px;
  color:#0EA5A5;
}

.workflow-arrow{
  color:#94a3b8;
  font-size:14px;
}
.landing-card{
  min-height:260px;
}
.content-wrapper,
.right-side,
.content {

  background:#ffffff !important;

}
.skin-blue .wrapper,
.skin-blue .main-sidebar,
.skin-blue .left-side {

  background-color:#ffffff;

}
.landing-description{
  font-size:1px;
  color:#64748b;
  max-width:700px;
  margin:0 auto;
  line-height:1.4;
}
.landing-card{

  height:260px;

  display:flex;

  flex-direction:column;

  background:white;

  border-radius:12px;

  padding:10px;

  box-shadow:0 2px 6px rgba(0,0,0,.08);

}
.landing-card img{

  height:145px;

  width:100%;

  object-fit:cover;

  border-radius:8px;

}
.landing-card h4{

  font-size:15px;

  margin-top:8px;

  margin-bottom:4px;

}
.landing-card p{

  font-size:12px;

  line-height:1.25;

  color:#64748b;

  margin-bottom:0;

}

.landing-page{

  height:calc(100vh - 50px);

  display:flex;

  flex-direction:column;

  justify-content:space-evenly;

  overflow:hidden;

}
.landing-brand{

  display:flex;

  align-items:center;

  justify-content:center;

  gap:16px;

}
/* ==========================
   Font
   ========================== */
body,
.main-header,
.main-sidebar,
.content-wrapper,
.box,
.form-control,
.btn,
.dataTable,
table {

    font-family: var(--font-family) !important;

}


    h1, h2, h3, h4, h5 {
    font-weight: 600;
    }



    .sidebar-menu > li > a {
        font-weight: 500;
    }

    .box-body {
    line-height: 1.5;
      }

      .control-label {
          font-weight: 500;
      }

/* ==========================
   SIDEBAR
   ========================== */

.main-sidebar {
    background-color: var(--sidebar-bg) !important;
}

.sidebar-menu > li > a {
    color: var(--text-light) !important;
}

.sidebar-menu > li:hover > a,
.sidebar-menu > li.active > a {

    background-color: rgba(255,255,255,0.06) !important;

    color: white !important;

    border-left: 3px solid var(--color-accent) !important;

}

#sidebarHandle {

    color: white;

}
#sidebarHandle {

    width: 100%;

    background: transparent;

    padding-top: 8px !important;

    padding-bottom: 8px !important;

    transition: all 0.2s ease;

}
#sidebarHandle i {

    color: var(--sidebar-icon) !important;

}


#sidebarHandle:hover i {

    color: white !important;

}

.sidebar-collapse #sidebarHandle {

    text-align: center !important;

    padding-left: 0 !important;

    padding-right: 0 !important;

}
.sidebar-menu {

    margin-top: 0 !important;

}

/* ==========================
   HEADER
   ========================== */

.main-header .logo {

    width: 360px !important;


    width: 100% !important;

    background-color: var(--header-bg) !important;

    display: flex !important;

    align-items: center !important;

justify-content: flex-start !important;

padding-left: 20px !important;


    overflow: hidden;
    height: 50px !important;


}

.main-header .navbar {

    display: none !important;

}




/* ==========================
   Menu spacing
   ========================== */
.sidebar-menu > li > a {

    padding-top: 14px !important;
    padding-bottom: 14px !important;

    font-size: 14px;

    border-radius: 8px;

    margin-left: 8px;
    margin-right: 8px;

}

/* ==========================
   Icons. 
   ========================== */
.sidebar-menu .fa,
.sidebar-menu .fas,
.sidebar-menu .far,
.sidebar-menu .fab {

    width: 22px;

    text-align: center;

    margin-right: 10px;

    font-size: 16px;

    color: var(--sidebar-icon);

    transition: all 0.2s ease;

}
.sidebar-menu > li.active i {

    color: var(--color-accent) !important;

}
.sidebar-menu > li:hover i {

    color: white !important;

}

/* ==========================
   Fondo
   ========================== */
.content-wrapper,
.right-side {

    background-color: var(--page-bg) !important;

}


/* ==========================
   boxes
   ========================== */

.box {

    background: var(--card-bg) !important;

    border-radius: 16px !important;

    overflow: visible !important;

    border: none !important;

box-shadow:
    0 4px 12px rgba(28,16,68,0.06),
    0 1px 3px rgba(0,0,0,0.04) !important;

    transition:
        transform 0.2s ease,
        box-shadow 0.2s ease;

}
.box-header {

    border-top-left-radius: 16px !important;
    border-top-right-radius: 16px !important;

}

.box-body {

    border-bottom-left-radius: 16px !important;
    border-bottom-right-radius: 16px !important;

}

.box:hover {

    transform: translateY(-2px);

    box-shadow:
        0 8px 20px rgba(0,0,0,0.08) !important;

}

/* All card types look identical */

.box-primary,
.box-info,
.box-warning,
.box-success {

    border: none !important;

    border-top: 2px solid var(--color-secondary) !important;

}

/* Header */

.box-header {

    background: var(--card-header-bg) !important;

    color: var(--text-primary) !important;

    border-bottom: 1px solid var(--card-border) !important;

    padding-top: 14px !important;

    padding-bottom: 14px !important;

}


/* Title */

.box-title {

    font-size: 15px !important;

    font-weight: 600 !important;

    color: var(--text-primary) !important;

}

/* Collapse / expand icons */

.btn-box-tool {

    color: var(--text-secondary) !important;

    font-size: 14px !important;

}

.btn-box-tool:hover {

    color: var(--color-accent) !important;

}

/* ==========================
   QC INFO CARDS
   ========================== */

.qc-info-card {

    background: var(--card-bg);

    border-radius: 16px;

    padding: 22px;

    margin-bottom: 10px;

    border-top: none;

    box-shadow:
        0 4px 12px rgba(28,16,68,0.06),
        0 1px 3px rgba(0,0,0,0.04);

}

.qc-info-title {

    font-size: 12px;

    text-transform: uppercase;

    letter-spacing: 0.08em;

    color: var(--text-secondary);

    margin-bottom: 8px;

}

.qc-info-value {

    font-size: 34px;

    font-weight: 700;

    color: var(--text-primary);

    line-height: 1;

}

.qc-info-subtitle {

    margin-top: 8px;

    color: var(--text-secondary);

    font-size: 13px;

}
/* ==========================
   controles
   ========================== */
.form-control {

    border-radius: 10px !important;

    border: 1px solid var(--input-border) !important;

    box-shadow: none !important;

    transition: all 0.2s ease;

}
.form-group {

    margin-bottom: 18px !important;

}

.form-control:focus {

    border-color: var(--color-primary) !important;

    box-shadow: 0 0 0 3px var(--color-primary-ring) !important;

}
.control-label {

    font-size: 13px;

    text-transform: uppercase;

    letter-spacing: 0.04em;

    color: var(--text-secondary);

}
input[type='checkbox']:checked {

    accent-color: var(--color-primary);

}
.selectize-input {

    border-radius: 10px !important;

    border: 1px solid var(--input-border) !important;

}

/* Radio buttons */

input[type='radio'] {

    accent-color: var(--color-primary);

}

/* ==========================
  Botones
   ========================== */

.btn {

    border-radius: 10px !important;

    font-weight: 500 !important;

}

.btn-default {

    background-color: var(--color-accent) !important;

    border-color: var(--color-accent) !important;

    color: white !important;

}

.btn-default:hover {

    background-color: var(--color-accent-hover) !important;

    border-color: var(--color-accent-hover) !important;

    color: white !important;

}

/* ==========================
  tablas
   ========================== */



   table.dataTable {

    border-collapse: separate !important;

    border-spacing: 0 6px !important;

}
table.dataTable thead th {

    background-color: var(--table-header) !important;

    border-bottom: 1px solid var(--table-border) !important;


    font-weight: 600 !important;

}
table.dataTable thead th {

    color: white !important;

}
table.dataTable tbody tr {

    background-color: white;

}

table.dataTable tbody tr:hover,
table.dataTable tbody tr:hover > td,
table.dataTable.display tbody tr:hover > .sorting_1,
table.dataTable.hover tbody tr:hover > .sorting_1 {

    background-color: var(--table-hover) !important;
    color: var(--text-primary) !important;

}

table.dataTable,
table.dataTable th,
table.dataTable td {

    font-family: 'Segoe UI', sans-serif !important;

}
.dataTables_wrapper {

    font-size: 13px;

}
table.dataTable tbody td {

    padding: 10px 12px !important;

}
table.dataTable tbody > tr.selected,
table.dataTable tbody > tr.selected > td,
table.dataTable tbody > tr > .selected,

table.dataTable tbody > tr.selected:hover,
table.dataTable tbody > tr.selected:hover > td,

table.dataTable.display tbody > tr.selected:hover > .sorting_1,
table.dataTable.hover tbody > tr.selected:hover > .sorting_1 {

    background-color: rgba(14,165,165,0.12) !important;
    color: var(--text-primary) !important;
    box-shadow: none !important;

}
table.dataTable tbody tr:hover:not(.selected) {

    background-color: var(--table-hover) !important;

}


/* ==========================
  help icons
   ========================== */
.fa-question-circle,
.fa-info-circle,
.fa-exclamation-circle {

    color: var(--text-secondary) !important;

    transition: color 0.2s ease;

}

.fa-question-circle:hover,
.fa-info-circle:hover,
.fa-exclamation-circle:hover {

  color: var(--color-accent) !important;

}
.nav-tabs > li.active .fa-question-circle,
.nav-tabs > li.active .fa-info-circle {

    color: rgba(255,255,255,0.8) !important;

}
.nav-tabs > li.active .fa-question-circle:hover,
.nav-tabs > li.active .fa-info-circle:hover {

    color: var(--color-accent) !important;

}


.tooltip-inner {

   background-color: var(--tooltip-bg) !important;

    color: white !important;

    border-radius: 10px !important;

    padding: 12px !important;

    font-size: 13px !important;

    line-height: 1.5 !important;

}
.tooltip {

    z-index: 999999 !important;

}
a.action-button:hover i,
a.action-button:hover .fa,
a.action-button:hover .fas {

    color: var(--color-accent) !important;

}
.modal {

    z-index: 500000 !important;

}

.modal-dialog {

    z-index: 500001 !important;

}

.modal-backdrop {

    z-index: 499999 !important;

}
.shiny-notification {

    z-index: 600000 !important;

}
/* ==================================
   Reference protein card
   ================================== */

.ref-protein-card {

    background: var(--card-bg);

    border: 1px solid var(--card-border);

    border-radius: 12px;

    padding: 14px;

    min-height: 100px;

}

.ref-protein-title {

    font-weight: 600;

    color: var(--text-primary);

    margin-bottom: 10px;

}

.ref-protein-title i {

    color: var(--color-accent);

    margin-right: 6px;

}

/* ==========================
  Tabs
   ========================== */
/* ==========================
   Segmented tabs
   ========================== */

.nav-tabs {

    border-bottom: none !important;

    display: flex;

    gap: 6px;

    margin-bottom: 12px;

}


.nav-tabs > li {

    margin-bottom: 0 !important;

}

.nav-tabs > li > a {

    border: none !important;

    border-radius: 8px !important;

    background: var(--tab-bg) !important;


    color: var(--tab-text) !important;

    font-weight: 500 !important;

    padding: 10px 18px !important;

    transition: all 0.2s ease;

}

.nav-tabs > li > a:hover {

    background: var(--tab-hover) !important;

    color: var(--text-primary) !important;

}

.nav-tabs > li.active > a,
.nav-tabs > li.active > a:hover,
.nav-tabs > li.active > a:focus {

    background: var(--color-secondary) !important;

    color: white !important;

}
.nav-tabs > li.active > a {

    box-shadow:
        0 4px 12px rgba(45,20,90,0.25);

}
/* Selectize */

.selectize-input.focus {

    border-color: var(--color-primary) !important;

    box-shadow: 0 0 0 3px var(--color-primary-ring) !important;

}

/* Selectize active option */

.selectize-dropdown .active {

    background-color: var(--color-primary-soft) !important;

    color: var(--text-primary) !important;

}




.box {

    position: relative;

}

.box:has(.selectize-dropdown),
.box:has(.selectize-input.focus) {

    z-index: 10000 !important;

}



.selectize-control.dropdown-active {

    z-index: 200000 !important;

}

.selectize-dropdown {

    z-index: 200001 !important;

}

.irs-bar,
.irs-bar-edge,
.irs-single {

    background: var(--color-primary) !important;

    border-color: var(--color-primary) !important;

}

.irs-slider {

    border-color: var(--color-primary) !important;

    background: var(--color-primary) !important;

}
.irs-from,
.irs-to,
.irs-single {

    background: var(--color-primary) !important;

    color: white !important;

}
.irs-from:after,
.irs-to:after,
.irs-single:after {

    border-top-color: var(--color-primary) !important;

}



.progress-bar {

    background-color: var(--color-accent) !important;

}



/* shinycssloaders */







    /* Fixed header */
    .main-header {
      position: fixed !important;
      top: 0;
      width: 100%;
      z-index: 99999 !important;
    }
    .main-header .logo {

    z-index: 100000 !important;

}

    /* Fixed sidebar */
    .main-sidebar {
      position: fixed !important;
      top: 50px !important;
      left: 0 !important;
      height: calc(100vh - 50px) !important;
      overflow-y: auto !important;
    }

    /* Prevent content from hiding under header */
    .content-wrapper,
    .right-side {
      padding-top: 50px !important;
    }
.tooltip-inner.curve-fit-content {
  max-width: 700px !important;
  width: 700px !important;
}

/* Sidebar collapse handle */

/* Reopen handle when sidebar is collapsed */

.sidebar-collapse .main-sidebar {

    transform: translateX(-230px);

    overflow: visible !important;

}


.sidebar-collapse #sidebarHandle {

    position: absolute;

    right: -24px;

    top: 0;

    width: 24px;

    height: 50px;

    background: var(--color-accent);

    border-radius: 0 8px 8px 0;

    display: flex !important;

    align-items: center;

    justify-content: center;

    border-bottom: none !important;

}
.main-sidebar {

    padding-top: 0 !important;

}




  "))
),

tags$script(HTML("

  // ✅ Fix DT + plots when box expands
$(document).on('expanded.boxwidget shown.bs.collapse', function(e) {

  setTimeout(function() {

    if ($.fn.dataTable) {
      $.fn.dataTable.tables({ visible: true, api: true }).columns.adjust();
    }

    window.dispatchEvent(new Event('resize'));
    window.dispatchEvent(new Event('resize'));
    window.dispatchEvent(new Event('resize'));

  }, 500);

});


  // ✅ Custom handler (called from Shiny server)
  Shiny.addCustomMessageHandler('triggerResize', function(message) {

    setTimeout(function() {

      window.dispatchEvent(new Event('resize'));
      window.dispatchEvent(new Event('resize'));

    }, 300);
  });


  // ✅ Fix DT when window resizes (fullscreen issue)
  $(window).on('resize', function() {
    clearTimeout(window.dtResizeTimeout);
    window.dtResizeTimeout = setTimeout(function() {

      if ($.fn.dataTable) {
        $.fn.dataTable.tables({ visible: true, api: true }).columns.adjust();
      }

    }, 200);
  });



$(document).on('mouseenter', '[data-toggle=\"tooltip\"]', function() {

  $(this).tooltip({
    html: true
  });

});

$(document).on('click', '[data-toggle=\"popover\"]', function() {

  $(this).popover({
    html: true,
    trigger: 'focus',
    container: 'body'
  });

});


"))
,

  useShinyjs(),

  use_busy_spinner(spin = "fading-circle"),  # Enable shinyjs
  tabItems(
    tabItem(
      tabName = "welcome",
      box(
        title = "Welcome to TPP Analyst",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        uiOutput("welcome_text")
      )
    ),

 tabItem(
        tabName = "data_input",
        fluidRow(

box(
  title = tagList(
    "Try the application",
tags$span(
  `data-toggle` = "tooltip",
  `data-placement` = "right",
   `data-html` = "true",
title = HTML(
  "<b>Example dataset</b><br><br>

  Staurosporine thermal proteome profiling dataset from the landmark study by Savitski et al. (Science, 2014).<br><br>

  DOI: 10.1126/science.1255784<br><br>

  Already log2 scaled and normalized to lowest temperature<br><br>
  
  Loading the example dataset automatically populates all required input files."
),
  tags$i(
    class = "fa fa-question-circle",
    style = "margin-left:8px; color:#337ab7; cursor:pointer;"
  )
)
  ),
    width = 12,
    status = "info",
    solidHeader = TRUE,

    p(
      "You can either upload your own files below or load the example dataset."
    ),

    actionButton(
      "load_example",
      "Load example dataset",
      icon = icon("database")
    )
  )

),

      # Invisible anchor to catch focus and stop scrolling
      div(id = "top_of_data_input", tabindex = "-1"),
      

        fluidRow(
          box(
            title = "Upload data",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            fluidRow(
              column(width = 6,
              fileInput(
                "file2",

                tagList(
                  "Experimental design ",
                    tags$span(
                      `data-toggle` = "tooltip",
                      `data-placement` = "right",
                      `data-html` = "true",

                      title = HTML(
                        "<b>Required columns</b><br><br>
                        • <b>File Name</b>: must match the measurement columns.<br>
                        • <b>Condition</b>: experimental condition.<br>
                        • <b>Temperature</b>: temperature value.<br>
                        • <b>Replicate</b>: replicate identifier."
                      ),

                      tags$i(
                        class = "fa fa-question-circle",
                        style = "margin-left:8px; color:#337ab7; cursor:pointer;"
                      )
                    )
                ),
                    accept = c(
                      ".csv",
                      ".tsv",
                      ".txt",
                      ".xlsx",
                      ".xls"
                    )
                  ),

                  uiOutput("protein_upload_message"),

                  uiOutput("protein_upload_ui"),
                  checkboxInput(
                  "already_log2",
                  "Click here if the dataset is already log2 transformed",
                  value = FALSE
                ),
                checkboxInput(
                  "Already_lowestnorm",
                  "Click here if the dataset is already normalized to the lowest temperature",
                  value = FALSE
                ),
                  uiOutput("custom_org_loader_ui")
              ),
              column(
                width = 6,

                uiOutput("condition_organism_mapping_ui"),
                uiOutput("custom_org_mapping_ui")
              )
            )
          )
        ),
        fluidRow(

          box(
            id = "data_preview_box_outer",
            title = "Data preview",
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE,
            collapsed = FALSE,

            div(class = "box-body",
                uiOutput("dataset_selector2"),
                div(
                  id = "data_preview_box",
                  style = "",   # adjust to your liking
                  withSpinner(DTOutput("head"),color = "#0EA5A5")
                )
            )
          )
          ,
box(
  id = "dataset_splitting_box",
  title = tagList(
    "Dataset splitting",
tags$span(
  `data-toggle` = "tooltip",
  `data-placement` = "right",
   `data-html` = "true",
title = HTML(
  "<b>Dataset splitting</b><br><br>

  Select how datasets should be organized for all downstream analyses.<br><br>

  <b>Condition Only</b><br>
  Creates one dataset for each experimental condition.<br><br>

  <b>Condition + Organism</b><br>
  Creates separate datasets for each condition-organism combination.<br><br>

  This option is recommended for experiments containing multiple organisms."
),
  tags$i(
    class = "fa fa-question-circle",
    style = "margin-left:8px; color:#337ab7; cursor:pointer;"
  )
)
  ),
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE,
            collapsed = FALSE,

            fluidRow(
            column(
              width = 4,
              uiOutput("split_mode_ui")
            ),
              column(
                width = 8,
                selectInput(
                  "selected_split",
                  "Select dataset:",
                  choices = character(0)
                )
              )
            ),

            hr(),

            
            div(
              style = "overflow-x: auto; width: 100%;",
              DTOutput("split_table")
            )

          )
        )
      ),
tabItem(
  tabName = "QC",

  # Conditional tabsetPanel to switch between Datasets
  uiOutput("qc_tabs_ui")
),

tabItem(
  tabName = "normalization",
  
  # Same idea as QC
  uiOutput("norm_tabs_ui")
),
tabItem(
  tabName = "clustering",

  uiOutput("clustering_tabs_ui")
),
        tabItem(
          tabName = "model_fitting",
          uiOutput("model_fitting_tabs_ui")
      ),
 
tabItem(
  tabName = "protein_plots",
  fluidRow(
    box(
      title = "Protein selection",
      width = 12,
      status = "primary",
      solidHeader = TRUE,
      uiOutput("protein_plots_tabs_ui")   
    ),
    box(
      title = "Protein abundance across temperatures",
      width = 12,
      status = "info",
      solidHeader = TRUE,
      plotlyOutput("protein_plot", height = 600) %>% withSpinner(
color = "#0EA5A5"
)
    )
  )
),
tabItem(
  tabName = "go_terms",
  uiOutput("go_terms_tabs_ui")
),
tabItem(

  tabName = "condition_models",

  tabsetPanel(

    id = "condition_model_tabs",

    tabPanel(
      "NPARC",
      uiOutput("model_comparison_ui")
    ),

    tabPanel(
      "TPP",
      uiOutput("tpp_analysis_ui")
    )

  )
)
),



tags$head(
  tags$script(HTML("
    Shiny.addCustomMessageHandler('clearDTSelection', function(message) {
      var table = $('#'+message.id).find('table').DataTable();
      if (table) {
        table.rows().deselect();
      }
    });
  "))
)



)
      # JavaScript handler

ui <- dashboardPage(header, sidebar, body)