package view;

import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.fxml.FXML;
import javafx.geometry.Side;
import javafx.scene.Node;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Label;
import javafx.scene.control.Slider;
import umlKurs.ChartCreator;

import java.awt.*;

public class ChartController {

    public Slider timeSlider;
    public LineChart lineChart;
    public NumberAxis xAxis;
    public NumberAxis yAxis;
    public Label timeLabel;
    private ChartCreator chartCreator;
    private double[][] chartSeparArray;
    private double[][] chartRunArray;
    private XYChart.Series seriesSepar = new XYChart.Series();
    private XYChart.Series seriesRun = new XYChart.Series();
    private ObservableList<XYChart.Data<Number, Number>> dataSepar = FXCollections.observableArrayList();
    private ObservableList<XYChart.Data<Number, Number>> dataRun = FXCollections.observableArrayList();

    @FXML
    private void initialize() {
        lineChart.setLegendVisible(true);
        lineChart.setLegendSide(Side.TOP);
        timeSlider.valueProperty().addListener(new ChangeListener<Number>() {
            @Override public void changed(ObservableValue<? extends Number> observableValue, Number oldValue, Number newValue) {
                int i = (int) Math.round((double)newValue*chartCreator.getNt()/chartCreator.getTime());
                timeLabel.setText("t = " + i*chartCreator.getTime()/chartCreator.getNt());
                for(int j = 0; j <= chartCreator.getNx(); j++) {
                    double separValue = chartSeparArray[j][i];
                    double runValue = chartRunArray[j][i];
                    dataSepar.get(j).setYValue(separValue);
                    dataRun.get(j).setYValue(runValue);
                }
            }
        });
    }

    public void setChartCreator(ChartCreator chartCreator) {
        this.chartCreator = chartCreator;
        chartSeparArray = chartCreator.composeSeparChartArray();
        chartRunArray = chartCreator.composeRunChartArray();
        timeSlider.setMin(0);
        timeSlider.setMax(chartCreator.getTime());
        timeSlider.setValue(0);
        yAxis.setUpperBound(1);
        yAxis.setLowerBound(-0.15);
        yAxis.setAutoRanging(false);
        xAxis.setUpperBound(chartCreator.getLength()+5);
        xAxis.setLowerBound(-5);
        xAxis.setAutoRanging(false);
        lineChart.setCreateSymbols(false);
        lineChart.setLegendVisible(false);
        lineChart.getStyleClass().add("thick-chart");

        createChart();
    }

    private void createChart(){
        double length = chartCreator.getLength();
        double xStep = length/chartCreator.getNx();
        for(int j = 0; j <= chartCreator.getNx(); j++) {
            double separValue = chartSeparArray[j][0];
            double runValue = chartRunArray[j][0];
            double x = j*xStep;
            dataSepar.add(new XYChart.Data<Number, Number>(x, separValue));
            dataRun.add(new XYChart.Data<Number, Number>(x, runValue));
        }
        seriesSepar.setName("Separ");
        seriesRun.setName("Run");

        seriesSepar.getData().addAll(dataSepar);
        seriesRun.getData().addAll(dataRun);
        lineChart.getData().addAll(seriesSepar);
        lineChart.getData().addAll(seriesRun);

        Node lineSepar = seriesSepar.getNode().lookup(".chart-series-line");
        Color color = Color.RED;
        String rgb = String.format("%d, %d, %d",
                (color.getRed() * 255) - 30,
                (color.getGreen() * 255) - 30,
                (color.getBlue() * 255) - 30);
        lineSepar.setStyle("-fx-stroke: rgba(" + rgb + ", 1.0);");


        Node lineRun = seriesRun.getNode().lookup(".chart-series-line");
        color = Color.BLACK;
        rgb = String.format("%d, %d, %d",
                (color.getRed() * 255),
                (color.getGreen() * 255),
                (color.getBlue() * 255));
        lineRun.setStyle("-fx-stroke: rgba(" + rgb + ", 1.0);");
        timeLabel.setText("t = 0");
    }
}
