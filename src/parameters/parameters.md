---
title: ""
header-includes: \usepackage{caption}
date: ""
output:
  html_document:
    keep_md: yes
    self_contained: yes
    css: styles.css
    theme: cosmo
    fig_caption: TRUE
  pdf_document: 
    fig_caption: true
    includes: 
      before_body: footer.tex
  word_document:
    df_print: kable
    fig_caption: yes
    fig_height: 4
    fig_width: 6
    highlight: tango
    reference_docx: template.docx
---

<script type="text/javascript">

if ($(window).width() < 768) {
  $('.dropdown-menu a.dropdown-toggle').on('click', function(e) {
    if (!$(this).next().hasClass('show')) {
      $(this).parents('.dropdown-menu').first().find('.show').removeClass("show");
    }

 

    var $subMenu = $(this).next(".dropdown-menu");
    if (!$subMenu.hasClass('show')) {
      $subMenu.addClass('show');
      $subMenu.show();
    } else {
      $subMenu.removeClass('show');
      $subMenu.hide();
    }

 

    $(this).parents('li.nav-item.dropdown.show').on('hidden.bs.dropdown', function(e) {
      $('.dropdown-submenu .show').removeClass("show");
    });

 

    return false;
  });
}

</script>



## Parameters

<br> 
The parameter table below summarises the current best estimates incorporated in the package. These will be updated as our understanding of the epidemic develops.

| **Parameter** | **Value** | **Reference** |
| --- | --- | --- |
| Basic reproductive number, R0 | 3.0 | [Report 13](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-13-europe-npi-impact/) |
| Mean Incubation Period | 3 days | [Report 9](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-9-impact-of-npis-on-covid-19/); [Linton et al.](https://www.medrxiv.org/content/medrxiv/early/2020/01/28/2020.01.26.20018754.full.pdf); [Li et al.](https://www.nejm.org/doi/full/10.1056/NEJMoa2001316) |
| Generation Time | 6.75 days | [Report 9](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-9-impact-of-npis-on-covid-19/) |
| Mean Duration in I\_MILD | 2.1 days | Incorporates 0.5 days of infectiousness prior to symptoms; with parameters below ~95% of all infections are mild. In combination with mean duration in I\_CASE this gives a mean generation time as above |
| Mean Duration in I\_CASE | 4.5 days | Mean onset-to-admission of 4 days from UK data. Includes 0.5 days of infectiousness prior to symptom onset |
| Mean Duration of Hospitalisation for non-critical Cases (I\_HOSP) | 5 days | Based on unpublished UK data |
| Mean duration of Critical Care (I\_ICU) if survive | 7.3 days | Based on [UK data](https://www.icnarc.org/Our-Audit/Audits/Cmp/Reports) adjusted for censoring |
| Mean duration of Critical Care (I\_ICU) if die | 6 days | Based on [UK data](https://www.icnarc.org/Our-Audit/Audits/Cmp/Reports) |
| Mean duration of Stepdown post ICU (I\_Rec) | 2 days | Based on unpublished UK data |
| Mean duration of hospitalisation if require ICU but do not receive it | 1 day | Working assumption |
| Probability of dying in critical care | 50% | Based on [UK data](https://www.icnarc.org/Our-Audit/Audits/Cmp/Reports) |
| Probability of death if require critical care but do not receive it | 95% | Working assumption |
| Multiplier of hazard of death if require non-critical hospitalisation but oxygen/hospital bed is not available | 2 | Working assumption; doubles age-dependent hazard of death |

Age-Specific Parameters

| **Age-Group** | **Proportion of Infections Hospitalised** | **Proportion of hospitalised cases requiring critical care** | **Proportion of non-critical care cases dying** |
| --- | --- | --- | --- |
| 0 to 4 | 0.001 | 0.050 | 0.013 |
| 5 to 9 | 0.001 | 0.050 | 0.013 |
| 10 to 14 | 0.001 | 0.050 | 0.013 |
| 15 to 19 | 0.002 | 0.050 | 0.013 |
| 20 to 24 | 0.005 | 0.050 | 0.013 |
| 25 to 29 | 0.010 | 0.050 | 0.013 |
| 30 to 34 | 0.016 | 0.050 | 0.013 |
| 35 to 39 | 0.023 | 0.053 | 0.013 |
| 40 to 44 | 0.029 | 0.060 | 0.015 |
| 45 to 49 | 0.039 | 0.075 | 0.019 |
| 50 to 54 | 0.058 | 0.104 | 0.027 |
| 55 to 59 | 0.072 | 0.149 | 0.042 |
| 60 to 64 | 0.102 | 0.224 | 0.069 |
| 65 to 69 | 0.117 | 0.307 | 0.105 |
| 70 to 74 | 0.146 | 0.386 | 0.149 |
| 75 to 79 | 0.177 | 0.461 | 0.203 |
| 80+ | 0.180 | 0.709 | 0.580 |
| Source | [Verity et al. 2020](https://www.thelancet.com/pdfs/journals/laninf/PIIS1473-3099(20)30243-7.pdf) corrected for non-uniform attack rate in China (see [Report 12](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-12-global-impact-covid-19/)) | Adjusted from IFR distributional shape in [Verity et al. 2020](https://www.thelancet.com/pdfs/journals/laninf/PIIS1473-3099(20)30243-7.pdf) to give an overall proportion of cases requiring critical care of ~30% to match [UK data](https://www.icnarc.org/Our-Audit/Audits/Cmp/Reports) | Calculated from IFR in [Verity et al. 2020](https://www.thelancet.com/pdfs/journals/laninf/PIIS1473-3099(20)30243-7.pdf) corrected for non-uniform attack rate in China (see [Report 12](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-12-global-impact-covid-19/)) given the 50% fatality rate in critical care. |
<br>
