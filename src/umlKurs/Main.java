package umlKurs;

import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.layout.AnchorPane;
import javafx.stage.Stage;
import view.ChartController;
import view.MainViewController;

import java.io.IOException;
import java.sql.SQLException;


public class Main extends Application{
    private ChartController chartController;

    public void setChartCreator(ChartCreator chartCreator) {
        chartController.setChartCreator(chartCreator);
    }

    public void initChartForm() throws IOException, SQLException {
        FXMLLoader loader = new FXMLLoader(getClass().getResource("/view/Chart.fxml"));
        AnchorPane pane = loader.load();
        chartController = loader.getController();
        Stage stage = new Stage();
        stage.setTitle("График");
        stage.setScene(new Scene(pane));
        stage.show();
    }

    public static void main(String[] args) {
        launch(args);
    }

    @Override
    public void start(Stage primaryStage) throws Exception {
        FXMLLoader loader = new FXMLLoader(getClass().getResource("/view/MainView.fxml"));
        Parent root = loader.load();

        MainViewController startController = loader.getController();
        startController.setMain(this);
        primaryStage.setTitle("Параметры");
        primaryStage.setScene(new Scene(root));
        primaryStage.show();
    }
}
