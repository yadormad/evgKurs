package view;

import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.scene.control.TextField;
import umlKurs.ChartCreator;
import umlKurs.Main;

import java.io.IOException;
import java.sql.SQLException;

public class MainViewController {

    public TextField l;
    public TextField D;
    public TextField T;
    public TextField N;
    public TextField xStep;
    public TextField tStep;
    private Main main;

    @FXML
    private void initialize() { }

    public void setMain(Main main) {
        this.main = main;
    }

    public void buildChart(ActionEvent actionEvent) throws IOException, SQLException {
        ChartCreator chartCreator = new ChartCreator(Double.parseDouble(l.getText()),
                Double.parseDouble(D.getText()), Double.parseDouble(T.getText()), Integer.parseInt(N.getText()), Double.parseDouble(xStep.getText()), Double.parseDouble(tStep.getText()));
        main.initChartForm();
        main.setChartCreator(chartCreator);
    }
}
