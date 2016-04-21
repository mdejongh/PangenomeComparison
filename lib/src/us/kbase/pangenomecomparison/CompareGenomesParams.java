
package us.kbase.pangenomecomparison;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: CompareGenomesParams</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "pangenome_id",
    "pangenome_ws",
    "protcomp_id",
    "protcomp_ws",
    "output_id",
    "workspace"
})
public class CompareGenomesParams {

    @JsonProperty("pangenome_id")
    private String pangenomeId;
    @JsonProperty("pangenome_ws")
    private String pangenomeWs;
    @JsonProperty("protcomp_id")
    private String protcompId;
    @JsonProperty("protcomp_ws")
    private String protcompWs;
    @JsonProperty("output_id")
    private String outputId;
    @JsonProperty("workspace")
    private String workspace;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("pangenome_id")
    public String getPangenomeId() {
        return pangenomeId;
    }

    @JsonProperty("pangenome_id")
    public void setPangenomeId(String pangenomeId) {
        this.pangenomeId = pangenomeId;
    }

    public CompareGenomesParams withPangenomeId(String pangenomeId) {
        this.pangenomeId = pangenomeId;
        return this;
    }

    @JsonProperty("pangenome_ws")
    public String getPangenomeWs() {
        return pangenomeWs;
    }

    @JsonProperty("pangenome_ws")
    public void setPangenomeWs(String pangenomeWs) {
        this.pangenomeWs = pangenomeWs;
    }

    public CompareGenomesParams withPangenomeWs(String pangenomeWs) {
        this.pangenomeWs = pangenomeWs;
        return this;
    }

    @JsonProperty("protcomp_id")
    public String getProtcompId() {
        return protcompId;
    }

    @JsonProperty("protcomp_id")
    public void setProtcompId(String protcompId) {
        this.protcompId = protcompId;
    }

    public CompareGenomesParams withProtcompId(String protcompId) {
        this.protcompId = protcompId;
        return this;
    }

    @JsonProperty("protcomp_ws")
    public String getProtcompWs() {
        return protcompWs;
    }

    @JsonProperty("protcomp_ws")
    public void setProtcompWs(String protcompWs) {
        this.protcompWs = protcompWs;
    }

    public CompareGenomesParams withProtcompWs(String protcompWs) {
        this.protcompWs = protcompWs;
        return this;
    }

    @JsonProperty("output_id")
    public String getOutputId() {
        return outputId;
    }

    @JsonProperty("output_id")
    public void setOutputId(String outputId) {
        this.outputId = outputId;
    }

    public CompareGenomesParams withOutputId(String outputId) {
        this.outputId = outputId;
        return this;
    }

    @JsonProperty("workspace")
    public String getWorkspace() {
        return workspace;
    }

    @JsonProperty("workspace")
    public void setWorkspace(String workspace) {
        this.workspace = workspace;
    }

    public CompareGenomesParams withWorkspace(String workspace) {
        this.workspace = workspace;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((("CompareGenomesParams"+" [pangenomeId=")+ pangenomeId)+", pangenomeWs=")+ pangenomeWs)+", protcompId=")+ protcompId)+", protcompWs=")+ protcompWs)+", outputId=")+ outputId)+", workspace=")+ workspace)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
