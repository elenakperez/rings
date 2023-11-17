Is the Regime Shift in Gulf Stream Warm Core Rings Evident in Satellite Altimetry? An Inter-comparison of Eddy Identification and Tracking Products

Downstream of Cape Hatteras, the vigorously meandering Gulf Stream forms anticyclonic warm core rings (WCRs) that carry warm Gulf Stream and Sargasso Sea waters into the cooler, fresher Slope Sea, and forms cyclonic cold core rings (CCRs) that carry Slope Sea waters into the Sargasso Sea. The Northwest Atlantic shelf and open ocean off the U.S. East Coast have experienced dramatic changes in ocean circulation and water properties in recent years, with significant consequences for marine ecosystems and coastal communities. Some of these changes may be related to a reported regime shift in the number of WCRs formed annually, with a doubling of WCRs shed in the years after 2000. Since the regime shift was detected using a regional eddy-tracking product that is primarily based on sea surface temperatures from satellites and in-situ measurements and relies on analyst skill, we examine global eddy-tracking products as an automated and potentially more objective way to detect changes in Gulf Stream rings. Available global eddy tracking products rely on altimeter-measured sea surface height (SSH), with anticyclonic WCRs registering as sea surface highs and cyclonic CCRs as lows. To identify ocean eddies, these products use either contours of SSH or a Lagrangian approach with particles seeded in satellite-based surface geostrophic velocity fields (based on SSH gradients). Comparisons of these products with the regional SST-based product include annual WCR counts, seasonality, preferred formation location, regime shift detection, and several years of feature-to-feature WCR comparisons across products. This study confirms that the three global eddy products are not well suited for statistical analysis of rings shed from the Gulf Stream, a highly energetic region, and suggests that automated WCR identification and tracking comes at the price of accurate identification and tracking.

### Motivation: 
<a href="https://www.nature.com/articles/s41598-019-48661-9"> Gangopadhyay et al., 2019</a> found a regime shift in the number of warm core rings formations between 1980–1999 and 2000–2017. <a href="https://www.nature.com/articles/s41598-021-81827-y"> Silver et al., 2021</a> expanded upon this work and also found a regime shift in warm core rings, but *no* regime shift in cold core rings. **We seek to understand if a regime shift, and similar results, can be detected using a global mesoscale eddy dataset such as META2.0.** 

<div class="alert alert-block alert-info">
<b>The goal of this work is to:</b>
    <li>Repeat similar analyses as Silver et al., 2021 for eddy datasets that have been filtered to only include warm core ring- and cold core ring-like eddies</li>
    <li>Check if these eddy datasets detect a regime shift in Gulf Stream rings</li>
    <li>Analyize the formation patterns, seasonal cycles, and formation trends of the eddies</li>
    <li>Test if META2.0 and similar products could be used to study Gulf Stream rings</li>
</div>

<div class="alert alert-block alert-warning">
<b>Note:</b> The META2.0 dataset, or Chelton eddy tracks, are no longer recommended for use. Instead, AVISO reccomends using META3.1EXP. When this project began, that warning did not exist and META3.1EXP had not been published. Hence why we primarily use the META2.0 dataset.
</div>

### Datasets: 
☆<a href="https://www.aviso.altimetry.fr/en/data/products/value-added-products/global-mesoscale-eddy-trajectory-product/meta2-0-dt.html">The Altimetric Mesoscale Eddy Trajectories Atlas (META2.0) </a> <br>
The primary dataset used in this notebook is *The Altimetric Mesoscale Eddy Trajectories Atlas (META2.0)*. <a href="https://www.aviso.altimetry.fr/fileadmin/documents/data/products/value-added/Schlax_Chelton_2016.pdf"> Schlax & Chelton, 2016</a> give a product description of *META2.0*, which built upon <a href="https://www.sciencedirect.com/science/article/pii/S0079661111000036"> Chelton et al., 2011</a>. 

<a href="https://www.aviso.altimetry.fr/en/data/products/value-added-products/global-mesoscale-eddy-trajectory-product/meta3-1-exp-dt.html">The Mesoscale Eddy Trajectories Atlas (META3.1EXP) </a> <br>
The *Mesoscale Eddy Trajectories Atlas Product Verison 3.1 Experimental (META3.1EXP)* was first published in March 2022 and is an update to *META2.0*. The details of this product are described in <a href="https://essd.copernicus.org/articles/14/1087/2022/"> Pegliasco et al., 2022</a> which inherited an eddy-tracking algorithm developed by <a href="https://journals.ametsoc.org/view/journals/atot/31/5/jtech-d-14-00019_1.xml"> Mason et al., 2014</a> that is inpsired by <a href="https://www.sciencedirect.com/science/article/pii/S0079661111000036"> Chelton et al., 2011</a>.

<a href="https://zenodo.org/record/7349753">The Global Lagrangian Eddy Dataset (GLED v1.0) </a> <br>
The third eddy product differs from the *META* datasets because it uses Lagrangian methods to identify and track eddies. This dataset was created by <a href="https://essd.copernicus.org/articles/15/1765/2023/"> Liu and Abernathy, 2023</a> as the Lagrangian alternative to Eulerian eddy products like *META2.0* and all the successive Eulerian eddy products inspired by <a href="https://www.sciencedirect.com/science/article/pii/S0079661111000036"> Chelton et al., 2011</a>. 

<a href="https://www.bco-dmo.org/dataset/810182">Yearly census of Gulf Stream Warm Core Ring formation from 1980 to 2017</a><br>
The WCR census uses the <a href="https://jcgulfstream.com/charts/">Clark charts</a> to document the formation and demise dates, locations, and the area at formation for all WCRs formed between 1980 and 2017 that lived for a week or more. 
